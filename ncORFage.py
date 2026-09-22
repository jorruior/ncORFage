#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ncORFage.py
Compute Branch Length Scores (BLS) and ORF conservation ages
from a GTF of ORFs, using CodAlignView multi-species alignments and PRANK
ancestral sequence reconstruction.

Author: Jorge Ruiz-Orera
Date:   2026-04-01
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import time
import tempfile
from collections import defaultdict, namedtuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from io import StringIO
from math import floor
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
import warnings
from Bio import BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning, message="Partial codon")

try:
	from ete4 import Tree, PhyloTree
except ImportError:
	sys.exit("ERROR: ete4 is required.  Install with:  pip install ete4")
try:
	from Bio import Seq
except ImportError:
	sys.exit("ERROR: Biopython is required.  Install with:  pip install biopython")

try:
	import matplotlib
	matplotlib.use("Agg")
	matplotlib.rcParams["svg.fonttype"] = "none"
	import matplotlib.pyplot as plt
	from matplotlib.lines import Line2D
	from Bio import Phylo
	_PLOTTING_AVAILABLE = True
except ImportError:
	_PLOTTING_AVAILABLE = False

MAX_CODONS = 40000
CODALIGNVIEW_URL = "https://data.broadinstitute.org/compbio1/cav.php"
TREE_BASE_URL = "https://data.broadinstitute.org/compbio1/CodAlignViewFiles/TreeNHs/"
ALL_START_CODONS = {"ATG", "CTG", "GTG", "TTG"}
STOP_CODONS = {"TAG", "TGA", "TAA"}
LOG = logging.getLogger("orf_bls")

# Dot colors used in both tree plots
COLOR_CONSERVED = "#1f77b4"      # blue
COLOR_NOTCONSERVED = "#d62728"   # red
COLOR_MISSING = "#bbbbbb"        # grey (not aligned / invalid)



# Repeat overlap (RepeatMasker)

def parse_repeatmasker(rm_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
	"""
	Parse a RepeatMasker .out file into a dict: chrom -> sorted list of (start, end, class/family).
	Chromosome names are stored without 'chr' prefix for matching against GTF.
	"""
	repeats: Dict[str, list] = defaultdict(list)
	with open(rm_path) as fh:
		for line in fh:
			line = line.strip()
			if not line or line.startswith("SW") or line.startswith("score"):
				continue
			parts = line.split()
			if len(parts) < 11:
				continue
			try:
				chrom = parts[4].replace("chr", "")
				start = int(parts[5])
				end = int(parts[6])
				repeat_class = parts[10]  # class/family column
				repeats[chrom].append((start, end, repeat_class))
			except (ValueError, IndexError):
				continue
	# Sort by start position per chromosome
	for chrom in repeats:
		repeats[chrom].sort(key=lambda x: x[0])
	LOG.info("Loaded %d repeat elements from %s",
			 sum(len(v) for v in repeats.values()), rm_path)
	return dict(repeats)


def find_repeat_overlaps(
	chrom: str, exons: List[Tuple[int, int]],
	repeats: Dict[str, List[Tuple[int, int, str]]],
) -> str:
	"""
	Find repeat elements overlapping any exon of an ORF.
	Returns semicolon-separated list of repeat class/family names, or "NA".
	"""
	chrom_clean = chrom.replace("chr", "")
	rm_list = repeats.get(chrom_clean, [])
	if not rm_list:
		return "NA"

	import bisect
	starts = [r[0] for r in rm_list]
	max_ends = []
	for _, end, _ in rm_list:
		max_ends.append(max(end, max_ends[-1]) if max_ends else end)
	overlapping = set()

	for exon_start, exon_end in exons:
		# Search repeats starting before the exon ends; prefix maxima account for nested repeats.
		right = bisect.bisect_right(starts, exon_end)
		# Check backwards from there
		for i in range(right - 1, -1, -1):
			rs, re, rc = rm_list[i]
			if max_ends[i] < exon_start:
				break
			if rs <= exon_end and re >= exon_start:
				overlapping.add(rc)

	return ";".join(sorted(overlapping)) if overlapping else "NA"


# BLAST search

def setup_blast_db(fasta_path: str, outdir: str) -> str:
	"""
	Create a BLAST protein database from a FASTA file if not already present.
	Returns the path to the database prefix.
	"""
	db_prefix = os.path.join(outdir, "blastdb",
							 os.path.basename(fasta_path).replace(".fa", "").replace(".fasta", ""))
	os.makedirs(os.path.dirname(db_prefix), exist_ok=True)

	# Check if DB already exists
	if os.path.exists(db_prefix + ".pdb") or os.path.exists(db_prefix + ".psq"):
		return db_prefix

	cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot", "-out", db_prefix]
	LOG.info("Building BLAST database: %s", " ".join(cmd))
	result = subprocess.run(cmd, capture_output=True, text=True)
	if result.returncode != 0:
		raise RuntimeError(f"makeblastdb failed: {result.stderr[:500]}")
	return db_prefix


def run_blastp(query_seq: str, db_prefix: str, evalue: float, orf_id: str,
			   blast_dir: str) -> str:
	"""
	Run BLASTP for a single protein sequence against the database.
	Returns semicolon-separated list of hit subject IDs, "none" for no hits, or "failed" for search errors.
	"""
	# Write query to temp file
	query_file = os.path.join(blast_dir, f"{orf_id}_query.fa")
	with open(query_file, "w") as fh:
		fh.write(f">{orf_id}\n{query_seq}\n")

	out_file = os.path.join(blast_dir, f"{orf_id}_blast.tsv")

	cmd = [
		"blastp",
		"-query", query_file,
		"-db", db_prefix,
		"-evalue", str(evalue),
		"-outfmt", "6 sseqid evalue",
		"-max_target_seqs", "10",
		"-out", out_file,
	]
	try:
		result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
		if result.returncode != 0:
			LOG.debug("BLASTP failed for %s: %s", orf_id, result.stderr[:200])
			return "failed"

		if not os.path.exists(out_file):
			LOG.debug("BLASTP produced no output for %s", orf_id)
			return "failed"

		# Parse hits
		hits = set()
		if os.path.exists(out_file):
			with open(out_file) as fh:
				for line in fh:
					parts = line.strip().split("\t")
					if len(parts) >= 2:
						try:
							if float(parts[1]) <= evalue:
								hits.add(parts[0])
						except ValueError:
							continue
		# Clean up temp files
		for f in (query_file, out_file):
			try:
				os.remove(f)
			except OSError:
				pass

		return ";".join(sorted(hits)) if hits else "none"

	except (subprocess.TimeoutExpired, Exception) as exc:
		LOG.debug("BLASTP error for %s: %s", orf_id, exc)
		return "failed"

def parse_gtf(gtf_path: str) -> Dict[str, dict]:
	orfs: Dict[str, dict] = {}
	feature_types = {"CDS", "ORF"}
	with open(gtf_path) as fh:
		for line in fh:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			parts = line.split("\t")
			if len(parts) < 9:
				continue
			chrom, source, feature, start, end, score, strand, frame, attrs = parts[:9]
			if feature not in feature_types:
				continue
			orf_id = _parse_gtf_attr(attrs, "orf_id")
			if orf_id is None:
				orf_id = _parse_gtf_attr(attrs, "gene_id")
			if orf_id is None:
				continue
			start, end = int(start), int(end)
			if orf_id not in orfs:
				orfs[orf_id] = dict(chrom=chrom, strand=strand, exons=[])
			orfs[orf_id]["exons"].append((start, end))
	for orf in orfs.values():
		orf["exons"].sort(key=lambda x: x[0])
	LOG.info("Parsed %d ORFs from GTF %s", len(orfs), gtf_path)
	return orfs

def _parse_gtf_attr(attrs_str: str, key: str) -> Optional[str]:
	m = re.search(rf'{key}\s+"([^"]+)"', attrs_str)
	if m:
		return m.group(1)
	m = re.search(rf'{key}=([^;]+)', attrs_str)
	return m.group(1) if m else None



# CodAlignView alignment download

def _get_url(url: str, timeout: int = 60) -> str:
	try:
		with urlopen(url, timeout=timeout) as resp:
			return resp.read().decode("utf-8")
	except (HTTPError, URLError) as exc:
		raise RuntimeError(f"URL fetch failed for {url}: {exc}")
	except TimeoutError:
		raise RuntimeError(f"Timeout after {timeout}s for {url}")

def gtf_to_intervals(orf: dict) -> str:
	chrom = orf["chrom"]
	if not chrom.startswith("chr"):
		chrom = "chr" + chrom
	return "+".join(f"{chrom}:{s}-{e}" for s, e in orf["exons"])

def download_alignment(intervals: str, strand: str, alnset: str) -> str:
	strand_char = strand if strand in ("+", "-") else "+"
	url = (f"{CODALIGNVIEW_URL}?i={intervals}&s={strand_char}"
		   f"&a={alnset}&m={MAX_CODONS}&fo=")
	max_retries = 3
	last_exception = None
	for attempt in range(max_retries):
		try:
			fasta_str = _get_url(url, timeout=60)
			if not fasta_str.strip().startswith(">"):
				msg = re.sub(r"<[^>]*>", "", fasta_str).replace("CodAlignView Error\n", "")
				raise RuntimeError(f"CodAlignView error for {intervals}: {msg.strip()}")
			return fasta_str
		except (HTTPError, URLError, RuntimeError) as exc:
			last_exception = exc
			if attempt < max_retries - 1:
				LOG.debug("Attempt %d failed: %s. Waiting 30s...", attempt + 1, exc)
				time.sleep(30)
	raise RuntimeError(
		f"Could not download alignment for {intervals} {strand} after {max_retries} attempts: {last_exception}")

def download_alnset_tree(alnset: str) -> str:
	url = TREE_BASE_URL + alnset + ".nh"
	try:
		return _get_url(url, timeout=30).rstrip()
	except (HTTPError, URLError, RuntimeError):
		raise RuntimeError(f"Could not download tree for alnset '{alnset}'.")

def parse_fasta_str(fasta_str: str) -> List[Tuple[str, str]]:
	pairs = []
	name = None
	seq_parts: List[str] = []
	for line in fasta_str.split("\n"):
		line = line.rstrip()
		if not line:
			continue
		if line.startswith(">"):
			if name is not None:
				pairs.append((name, "".join(seq_parts)))
			name = line[1:].strip()
			seq_parts = []
		else:
			seq_parts.append(line)
	if name is not None:
		pairs.append((name, "".join(seq_parts)))

	# Keep the first sequence for each name; downstream lookups require unique names.
	seen = set()
	deduped = []
	for nm, sq in pairs:
		if nm in seen:
			continue
		seen.add(nm)
		deduped.append((nm, sq))
	return deduped

def _remove_allgap_columns(names, seqs):
	if not seqs:
		return list(zip(names, seqs))
	aln_len = len(seqs[0])
	gap_cols = set()
	for col in range(aln_len):
		if all(s[col] == '-' for s in seqs):
			gap_cols.add(col)
	if not gap_cols:
		return list(zip(names, seqs))
	return [(name, "".join(ch for i, ch in enumerate(seq) if i not in gap_cols))
			for name, seq in zip(names, seqs)]

def prepare_codon_alignment(fasta_str, ref_species, invalid_cutoff, min_species=2):
	pairs = parse_fasta_str(fasta_str)
	if not pairs:
		raise RuntimeError("Empty alignment")
	if not any(name == ref_species for name, _ in pairs):
		raise RuntimeError(f"Reference species '{ref_species}' not in alignment")

	iupac_nt = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
	precleaned = []
	for name, seq in pairs:
		if "." in seq or "|" in seq:
			continue
		seq = seq.replace("X", "N").replace("x", "n")
		seq = "".join(ch if (ch in iupac_nt or ch == "-") else "N" for ch in seq)
		ungapped = seq.replace("-", "")
		if len(ungapped) < 6:
			continue
		precleaned.append((name, seq))

	if not precleaned:
		raise RuntimeError("Reference species sequence too short after processing")

	ref_seq_full = next((s for n, s in precleaned if n == ref_species), None)
	if ref_seq_full is None:
		raise RuntimeError("Reference species sequence too short after processing")
	ref_seq_full_u = ref_seq_full.upper()

	filtered_names, filtered_seqs = [], []
	for name, seq in precleaned:
		ungapped = seq.replace("-", "")
		if len(ungapped) < 6:
			continue
		if name == ref_species:
			# Check reference N content without counting alignment gaps.
			if ungapped.upper().count("N") / len(ungapped) > invalid_cutoff:
				continue
		else:
			# Count gaps and Ns as invalid, excluding columns gapped in both species.
			seq_u = seq.upper()
			considered = 0
			n_invalid = 0
			for ref_ch, ch in zip(ref_seq_full_u, seq_u):
				if ref_ch == '-' and ch == '-':
					continue
				considered += 1
				if ch == '-' or ch == 'N':
					n_invalid += 1
			if considered == 0 or n_invalid / considered > invalid_cutoff:
				continue
		filtered_names.append(name)
		filtered_seqs.append(seq)

	if not any(name == ref_species for name in filtered_names):
		raise RuntimeError("Reference species sequence too short after processing")

	enough_species = len(filtered_names) >= min_species
	cleaned = _remove_allgap_columns(filtered_names, filtered_seqs)
	lines = []
	for name, seq in cleaned:
		lines.append(f">{name}")
		lines.append(seq)
	return "\n".join(lines) + "\n", enough_species



# PRANK ancestral reconstruction

def run_prank(fasta_path, tree_path, output_prefix):
	cmd = ["prank", f"-d={fasta_path}", f"-t={tree_path}", f"-o={output_prefix}",
		   "-showanc", "-keep", "-prunetree", "-DNA", "-once"]
	LOG.debug("Running: %s", " ".join(cmd))
	try:
		result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
	except subprocess.TimeoutExpired:
		raise RuntimeError(f"PRANK timed out after 300 seconds for {fasta_path}")
	if result.returncode != 0:
		raise RuntimeError(f"PRANK failed for {fasta_path}: {result.stderr[:500]}")
	asr_fas = output_prefix + ".anc.fas"
	asr_dnd = output_prefix + ".anc.dnd"
	if not os.path.exists(asr_fas):
		raise FileNotFoundError(f"PRANK did not produce expected file: {asr_fas}")
	if not os.path.exists(asr_dnd):
		raise FileNotFoundError(f"PRANK did not produce expected file: {asr_dnd}")
	_sed_inplace(asr_dnd, [(r'\)#', ')"Node'), (r'#', '"')])
	_sed_inplace(asr_fas, [(r'#', 'Node', 1), (r'#', '')])
	return asr_fas, asr_dnd

def _sed_inplace(filepath, replacements):
	with open(filepath, "r") as fh:
		lines = fh.readlines()
	new_lines = []
	for line in lines:
		for repl in replacements:
			pattern, sub = repl[0], repl[1]
			count = repl[2] if len(repl) > 2 else 0
			line = re.sub(pattern, sub, line, count=count)
		new_lines.append(line)
	with open(filepath, "w") as fh:
		fh.writelines(new_lines)



# ORF conservation logic

def remove_gaps(seq):
	return seq.replace("-", "").replace(".", "")

def to_codons(seq):
	return [seq[i:i+3] for i in range(0, len(seq), 3)]

def find_nth_nongap(seq, n, gap="-"):
	if n <= 0:
		return 0
	idx = count = 0
	while idx < len(seq) and count < n:
		if seq[idx] != gap:
			count += 1
		idx += 1
	return idx

ConservationResult = namedtuple(
	"ConservationResult", ["valid", "has_start", "no_pmsc", "conserved"])

def is_valid_sequence(seq, gapped_seq=None, invalid_cutoff=0.5, is_ref=False, ref_gapped_seq=None):
	seq = seq.upper()
	if '.' in seq:
		return False
	# Check invalid-content fraction against invalid_cutoff using the gapped
	# sequence when available, otherwise fall back to a check on ungapped.
	if gapped_seq is not None:
		gs = gapped_seq.upper()
		if '.' in gs:
			return False
		if is_ref:
			# For the reference, check N content in the ungapped sequence.
			if len(seq) > 0 and seq.count('N') / len(seq) > invalid_cutoff:
				return False
		elif ref_gapped_seq is not None and len(ref_gapped_seq) == len(gs):
			# Exclude shared gaps; count other gaps, Ns, and Xs as invalid.
			ref_gs = ref_gapped_seq.upper()
			considered = 0
			n_invalid = 0
			for ref_ch, ch in zip(ref_gs, gs):
				if ref_ch == '-' and ch == '-':
					continue
				considered += 1
				if ch in ('-', 'N', 'X'):
					n_invalid += 1
			if considered == 0 or n_invalid / considered > invalid_cutoff:
				return False
		else:
			# Count Ns, Xs, and gaps toward the invalid fraction.
			n_invalid = gs.count('N') + gs.count('-') + gs.count('X')
			if n_invalid / len(gs) > invalid_cutoff:
				return False
	else:
		if is_ref:
			if len(seq) > 0 and seq.count('N') / len(seq) > invalid_cutoff:
				return False
		else:
			n_invalid = seq.count('N') + seq.count('X')
			if len(seq) < 3:
				return len(seq) > 0 and n_invalid / len(seq) <= invalid_cutoff
			if n_invalid / len(seq) > invalid_cutoff:
				return False
	if len(seq) < 3:
		return len(seq) > 0
	if any(seq[i] == 'N' for i in (0, 1, 2, -1, -2, -3)):
		return False
	return True

def is_conserved(
	ref_seq, query_seq, stop_cutoff=0.7, invalid_cutoff=0.5,
	start_mode="same", initiation_offset=3, min_identity=0.0,
	query_is_ref=False,
):
	"""
	Check whether query_seq has conserved ORF structure relative to ref_seq.

	query_is_ref: set True when query_seq is itself the reference species'
	sequence (e.g. when computing the reference's own status for output) —
	in that case the query-side validity check is also exempted from gap/X
	filtering, same as the ref-side check.

	start_mode options:
	  "same"         : start codon must match the ref's first codon, or be ATG
	  "nearcognate"  : any of ATG/CTG/GTG/TTG within initiation_offset codons
	  "atg"          : only ATG within initiation_offset codons
	  "none"         : start codon not checked

	min_identity: minimum amino acid identity (0.0-1.0) vs ref to count as conserved.
	"""
	query_nogap = remove_gaps(query_seq).upper()
	ref_nogap = remove_gaps(ref_seq).upper()

	if not is_valid_sequence(ref_nogap, gapped_seq=ref_seq, invalid_cutoff=invalid_cutoff, is_ref=True):
		return ConservationResult(valid=False, has_start=False, no_pmsc=False, conserved=0)
	if not is_valid_sequence(query_nogap, gapped_seq=query_seq, invalid_cutoff=invalid_cutoff,
							 is_ref=query_is_ref, ref_gapped_seq=(None if query_is_ref else ref_seq)):
		return ConservationResult(valid=False, has_start=False, no_pmsc=False, conserved=0)

	query_codons = to_codons(query_nogap)
	ref_codons = to_codons(ref_nogap)

	# Start codon evaluation
	if start_mode == "none":
		has_start = True
		which_start = 0
	elif start_mode == "same":
		# Must match ref's first codon OR be ATG
		ref_first = ref_codons[0] if ref_codons else ""
		valid_starts = {ref_first, "ATG"} if ref_first else {"ATG"}
		check_codons = query_codons[:initiation_offset]
		hits = [c in valid_starts for c in check_codons]
		has_start = any(hits)
		which_start = hits.index(True) if has_start else 0
	else:
		valid_starts = ALL_START_CODONS if start_mode == "nearcognate" else {"ATG"}
		check_codons = query_codons[:initiation_offset]
		hits = [c in valid_starts for c in check_codons]
		has_start = any(hits)
		which_start = hits.index(True) if has_start else 0

	# Premature stop codon evaluation
	ref_aa_len = len(ref_nogap) // 3 - 1
	aa_required = floor(ref_aa_len * stop_cutoff)
	nt_required = (which_start + aa_required) * 3
	query_start_aln = find_nth_nongap(query_seq, which_start * 3)
	ref_end_aln = find_nth_nongap(ref_seq, nt_required)
	region = remove_gaps(query_seq[query_start_aln:ref_end_aln]).upper()

	if len(region) < 3:
		no_pmsc = False
	else:
		translated = Seq.Seq(region).translate()
		no_pmsc = "*" not in str(translated)

	# Identity check
	if min_identity > 0:
		aa_id = compute_aa_identity(ref_seq, query_seq)
		passes_identity = aa_id is not None and aa_id >= min_identity
	else:
		passes_identity = True

	# Combine
	if start_mode == "none":
		conserved = int(no_pmsc and passes_identity)
	else:
		conserved = int(has_start and no_pmsc and passes_identity)

	return ConservationResult(valid=True, has_start=has_start, no_pmsc=no_pmsc, conserved=conserved)


def compute_aa_identity(ref_seq, query_seq):
	ref_nogap = remove_gaps(ref_seq).upper()
	query_nogap = remove_gaps(query_seq).upper()
	if len(ref_nogap) < 3 or len(query_nogap) < 3:
		return None
	ref_aa = str(Seq.Seq(ref_nogap).translate())
	query_aa = str(Seq.Seq(query_nogap).translate())
	min_len = min(len(ref_aa), len(query_aa))
	if min_len == 0:
		return None
	total = sum(1 for i in range(min_len) if ref_aa[i] != '*')
	if total == 0:
		return None
	matches = sum(1 for i in range(min_len) if ref_aa[i] == query_aa[i] and ref_aa[i] != '*')
	return matches / total



# BLS calculation

def blsum(tree):
	return sum(n.dist for n in list(tree.traverse())[1:])

def find_most_distant_species(asr_tree_path, ref_sp, species_list):
	if not species_list or species_list == [ref_sp]:
		return ref_sp
	try:
		tree = Tree(open(asr_tree_path), parser=1)
		max_dist = -1
		most_distant = ref_sp
		for sp in species_list:
			if sp == ref_sp:
				continue
			try:
				d = tree.get_distance(ref_sp, sp)
				if d > max_dist:
					max_dist = d
					most_distant = sp
			except Exception:
				continue
		return most_distant
	except Exception:
		return "NA"

def _find_parent_node(tfull, leaf_name):
	"""Find the parent internal node of a leaf in the full tree."""
	try:
		leaf = tfull[leaf_name]
		parent = leaf.up
		if parent is not None:
			return parent.name
	except Exception:
		pass
	return "NA"

def _fixation_from_asr_node(tasr, origin_node_name, total_universe):
	"""
	Fraction of species/ancestors after the origin node in which the ORF is
	conserved (orf==1 and valid).

	total_universe: total number of descendants (leaves + internal nodes,
	excluding the origin node itself) in the full species tree subtree rooted
	at the equivalent global node.  Nodes absent from the ASR tree (unaligned
	species) count as non-conserved, so the denominator is max(tasr_total,
	total_universe).

	Returns a float rounded to 4 dp, or "NA" if the origin node is not found.
	"""
	try:
		origin_node = tasr[origin_node_name]
	except Exception:
		return "NA"
	conserved_total = 0
	asr_total = 0
	for n in origin_node.traverse():
		if n.name == origin_node_name:
			continue
		asr_total += 1
		if n.get_prop("valid") and n.get_prop("orf") == 1:
			conserved_total += 1
	# Unaligned species/ancestors are in tfull but absent from tasr; they count
	# as non-conserved, so use the larger of the two totals as denominator.
	denominator = max(asr_total, total_universe)
	if denominator == 0:
		return "NA"
	return round(conserved_total / denominator, 4)


def compute_bls(
	asr_fasta, asr_tree, full_tree_path, ref_sp,
	stop_cutoff=0.7, invalid_cutoff=0.5, start_mode="same",
	initiation_offset=3, min_identity=0.0,
):
	tasr = PhyloTree(open(asr_tree), parser=1, alignment=asr_fasta, alg_format="fasta")

	try:
		ref_node = tasr[ref_sp]
	except Exception:
		return _empty_result(ref_sp, reason="ref_species_not_in_tree")

	ref_seq = ref_node.props["sequence"]

	for n in tasr.traverse():
		cons = is_conserved(
			ref_seq, n.props["sequence"],
			stop_cutoff=stop_cutoff, invalid_cutoff=invalid_cutoff,
			start_mode=start_mode, initiation_offset=initiation_offset,
			min_identity=min_identity, query_is_ref=(n.name == ref_sp),
		)
		tasr[n.name].add_prop("orf", cons.conserved)
		tasr[n.name].add_prop("valid", cons.valid)

	if not tasr[ref_sp].get_prop("valid"):
		return _empty_result(ref_sp, reason="ref_invalid")

	# AA identity stats, and stop-codon stats among conserved sequences
	id_cons_vals, id_notcons_vals = [], []
	stop_cons_hits = []
	for n in tasr.leaves():
		if n.name == ref_sp or not n.get_prop("valid"):
			continue
		aa_id = compute_aa_identity(ref_seq, n.props["sequence"])
		if aa_id is None:
			continue
		if n.get_prop("orf") == 1:
			id_cons_vals.append(aa_id)
			n_nogap = remove_gaps(n.props["sequence"]).upper()
			stop_cons_hits.append(1 if n_nogap[-3:] in STOP_CODONS else 0)
		else:
			id_notcons_vals.append(aa_id)

	identity_conserved = round(sum(id_cons_vals)/len(id_cons_vals), 4) if id_cons_vals else "NA"
	identity_notconserved = round(sum(id_notcons_vals)/len(id_notcons_vals), 4) if id_notcons_vals else "NA"
	stop_conserved = round(sum(stop_cons_hits)/len(stop_cons_hits), 4) if stop_cons_hits else "NA"

	# Origination
	anc_names = [n.name for n in tasr[ref_sp].ancestors()]
	anc_orf = [n.get_prop("orf") for n in tasr[ref_sp].ancestors()]
	anc_valid = [n.get_prop("valid") for n in tasr[ref_sp].ancestors()]

	if not anc_names:
		r = _single_species_result(ref_sp)
		r["identity_conserved"] = identity_conserved
		r["identity_notconserved"] = identity_notconserved
		r["stop_conserved"] = stop_conserved
		return r

	if anc_orf[-1] == 0 or not anc_valid[-1]:
		origin_manner = "denovo"
		origin_age_local = ref_sp
		for i, (cons, val) in enumerate(zip(anc_orf[::-1], anc_valid[::-1])):
			if cons == 1 and val:
				origin_age_local = anc_names[-1 - i]
				break
	else:
		origin_manner = "nondenovo"
		origin_age_local = anc_names[-1]

	# Leaf sets (from ASR/local tree)
	origin_leaves = [n.name for n in tasr[origin_age_local].leaves() if n.get_prop("valid")]
	origin_leaves_orf_local = [
		n.name for n in tasr[origin_age_local] if n.get_prop("orf") == 1 and n.get_prop("valid")]
	local_leaves_orf = [
		n.name for n in tasr if n.get_prop("orf") == 1 and n.get_prop("valid")]

	if not origin_leaves_orf_local:
		origin_leaves_orf_local = [ref_sp]
	if not local_leaves_orf:
		local_leaves_orf = [ref_sp]

	# Full tree
	tfull = Tree(open(full_tree_path), parser=1)
	_name_internal_nodes(tfull)
	bl_all = blsum(tfull)

	# Normalize leaf names to match alignment and full-tree labels.
	def _norm_name(s):
		return "".join(c if c.isalnum() else "_" for c in s).lower()

	tfull_leaf_map = {_norm_name(n.name): n.name for n in tfull.leaves()}

	def _resolve(names):
		"""Map a list of ASR leaf names to their tfull equivalents."""
		out = []
		for nm in names:
			if nm in tfull_leaf_map.values():
				out.append(nm)
			else:
				mapped = tfull_leaf_map.get(_norm_name(nm))
				if mapped:
					out.append(mapped)
				# else: silently drop — species genuinely absent from full tree
		return out

	origin_leaves            = _resolve(origin_leaves)
	origin_leaves_orf_local  = _resolve(origin_leaves_orf_local)
	local_leaves_orf         = _resolve(local_leaves_orf)

	# Method 1: ASR-based origination node
	try:
		origin_age_global = tfull.common_ancestor(origin_leaves).name
	except Exception:
		r = _empty_result(ref_sp, reason="origin_leaves_not_in_full_tree")
		r["identity_conserved"] = identity_conserved
		r["identity_notconserved"] = identity_notconserved
		r["stop_conserved"] = stop_conserved
		return r

	bl_sub_origin = blsum(tfull[origin_age_global])

	if origin_leaves_orf_local == [ref_sp]:
		bl_orf_origin = tfull.get_distance(origin_age_global, ref_sp)
	else:
		try:
			mrca_name = tfull.common_ancestor(origin_leaves_orf_local).name
			t_pruned = tfull.copy()
			t_pruned.prune(origin_leaves_orf_local, preserve_branch_length=True)
			bl_orf_origin = tfull.get_distance(origin_age_global, mrca_name) + blsum(t_pruned)
		except Exception:
			bl_orf_origin = tfull.get_distance(origin_age_global, ref_sp)

	bls_global_origin = bl_orf_origin / bl_all if bl_all else 0
	bls_local_origin = bl_orf_origin / bl_sub_origin if bl_sub_origin else 0

	# Method 2: Naive MRCA — always a node in the full tree
	if local_leaves_orf == [ref_sp]:
		# Single species: use parent node as the naive age
		naive_age_global = _find_parent_node(tfull, ref_sp)
		bl_sub_naive = 0
		bl_orf_naive = 0
	else:
		try:
			naive_age_global = tfull.common_ancestor(local_leaves_orf).name
			bl_sub_naive = blsum(tfull[naive_age_global])
			t_naive = tfull[naive_age_global].copy()
			t_naive.prune(local_leaves_orf, preserve_branch_length=True)
			bl_orf_naive = blsum(t_naive)
		except Exception:
			naive_age_global = _find_parent_node(tfull, ref_sp)
			bl_sub_naive = 0
			bl_orf_naive = 0

	bls_global_naive = bl_orf_naive / bl_all if bl_all else 0
	bls_local_naive = bl_orf_naive / bl_sub_naive if bl_sub_naive else 0

	origin_species = sorted(n.name for n in tfull[origin_age_global].leaves())
	most_distant = find_most_distant_species(asr_tree, ref_sp, local_leaves_orf)

	# Fixation: fraction of species/ancestors descending from the origin node
	# in which the ORF is conserved.  Unaligned species (absent from tasr) count
	# as non-conserved; universe size is taken from the equivalent tfull subtree.

	# origin method
	origin_universe = sum(1 for n in tfull[origin_age_global].traverse()
	                      if n.name != origin_age_global)
	origin_age_global_fixation = _fixation_from_asr_node(
		tasr, origin_age_local, origin_universe)

	# naive method: MRCA of local_leaves_orf in tasr is the anchor.
	if local_leaves_orf == [ref_sp]:
		naive_age_global_fixation = "NA"
	else:
		try:
			naive_origin_local = tasr.common_ancestor(local_leaves_orf).name
			naive_universe = sum(1 for n in tfull[naive_age_global].traverse()
			                     if n.name != naive_age_global)
			naive_age_global_fixation = _fixation_from_asr_node(
				tasr, naive_origin_local, naive_universe)
		except Exception:
			naive_age_global_fixation = "NA"

	return {
		"origin_manner": origin_manner,
		"origin_age_local": origin_age_local,
		"origin_age_global": origin_age_global,
		"origin_species": ";".join(origin_species),
		"n_species_with_orf_local": len(local_leaves_orf),
		"species_with_orf_local": ";".join(sorted(local_leaves_orf)),
		"most_distant_orf_species_local": most_distant,
		"identity_conserved": identity_conserved,
		"identity_notconserved": identity_notconserved,
		"stop_conserved": stop_conserved,
		"bl_all": bl_all, "bl_sub_origin": bl_sub_origin,
		"bl_orf_origin": bl_orf_origin,
		"bls_global_origin": bls_global_origin, "bls_local_origin": bls_local_origin,
		"origin_age_global_fixation": origin_age_global_fixation,
		"naive_age_global": naive_age_global,
		"bl_sub_naive": bl_sub_naive, "bl_orf_naive": bl_orf_naive,
		"bls_global_naive": bls_global_naive, "bls_local_naive": bls_local_naive,
		"naive_age_global_fixation": naive_age_global_fixation,
	}

def _name_internal_nodes(tree):
	i = 1
	for node in tree.traverse("preorder"):
		if not node.is_leaf:
			node.name = f"N{i}"
			i += 1

def _empty_result(ref_sp, reason="unknown"):
	return {
		"origin_manner": reason, "origin_age_local": "NA", "origin_age_global": "NA",
		"origin_species": "NA", "n_species_with_orf_local": 0,
		"species_with_orf_local": ref_sp, "most_distant_orf_species_local": "NA",
		"identity_conserved": "NA", "identity_notconserved": "NA", "stop_conserved": "NA",
		"bl_all": "NA", "bl_sub_origin": "NA", "bl_orf_origin": "NA",
		"bls_global_origin": "NA", "bls_local_origin": "NA",
		"origin_age_global_fixation": "NA",
		"naive_age_global": "NA", "bl_sub_naive": "NA", "bl_orf_naive": "NA",
		"bls_global_naive": "NA", "bls_local_naive": "NA",
		"naive_age_global_fixation": "NA",
	}

def _single_species_result(ref_sp):
	return _empty_result(ref_sp, reason="single_species")

def _too_few_species_result(ref_sp):
	r = _empty_result(ref_sp, reason="too_few_species")
	r["origin_age_local"] = ref_sp
	r["origin_age_global"] = ref_sp
	r["origin_species"] = ref_sp
	r["n_species_with_orf_local"] = 1
	r["most_distant_orf_species_local"] = ref_sp
	return r

def _failed_result(ref_sp, reason="prank_failed"):
	r = _empty_result(ref_sp, reason=reason)
	r["origin_age_local"] = reason
	r["origin_age_global"] = reason
	return r



# Per-ORF alignment output

def compute_sequence_status(ref_seq, seq, stop_cutoff, invalid_cutoff, start_mode,
							initiation_offset, min_identity, is_ref=False):
	cons = is_conserved(ref_seq, seq, stop_cutoff=stop_cutoff,
						invalid_cutoff=invalid_cutoff, start_mode=start_mode,
						initiation_offset=initiation_offset, min_identity=min_identity,
						query_is_ref=is_ref)
	if not cons.valid:
		return -1, False
	return 1 if cons.conserved else 0, True

def translate_sequence(seq):
	ungapped = remove_gaps(seq).upper()
	if len(ungapped) < 3:
		return ""
	return str(Seq.Seq(ungapped).translate())

def write_orf_alignments(asr_fasta, asr_tree, orf_id, outdir, ref_sp,
						 stop_cutoff, invalid_cutoff, start_mode,
						 initiation_offset, min_identity):
	pairs = parse_fasta_str(open(asr_fasta).read())
	species_seqs = {name: seq for name, seq in pairs if not name.startswith("Node")}
	asr_seqs = {name: seq for name, seq in pairs if name.startswith("Node")}
	if not species_seqs:
		return
	ref_seq = species_seqs.get(ref_sp, "")
	if not ref_seq:
		return
	status = {}
	for name, seq in species_seqs.items():
		st, _ = compute_sequence_status(ref_seq, seq, stop_cutoff, invalid_cutoff,
										start_mode, initiation_offset, min_identity,
										is_ref=(name == ref_sp))
		status[name] = st
	asr_status = {}
	for name, seq in asr_seqs.items():
		st, _ = compute_sequence_status(ref_seq, seq, stop_cutoff, invalid_cutoff,
										start_mode, initiation_offset, min_identity)
		asr_status[name] = st

	# Build node -> descendant leaves map from the PRANK tree
	asr_node_leaves = {}
	if asr_seqs:
		try:
			tasr = Tree(open(asr_tree), parser=1)
			for node in tasr.traverse():
				if node.is_leaf:
					continue
				if node.name and node.name.startswith("Node"):
					leaves = sorted(l.name for l in node.leaves())
					asr_node_leaves[node.name] = ";".join(leaves)
		except Exception:
			pass

	with open(os.path.join(outdir, f"{orf_id}.aln.fa"), "w") as fh:
		for name, seq in species_seqs.items():
			fh.write(f">{name}|status={status[name]}\n{seq}\n")
	with open(os.path.join(outdir, f"{orf_id}.conserved.fa"), "w") as fh:
		fh.write(f">{ref_sp}|status={status[ref_sp]}\n{species_seqs[ref_sp]}\n")
		for name, seq in species_seqs.items():
			if name != ref_sp and status[name] == 1:
				fh.write(f">{name}|status=1\n{seq}\n")
	with open(os.path.join(outdir, f"{orf_id}.prot.fa"), "w") as fh:
		for name, seq in species_seqs.items():
			prot = translate_sequence(seq)
			if prot:
				fh.write(f">{name}|status={status[name]}\n{prot}\n")
	with open(os.path.join(outdir, f"{orf_id}.asr.fa"), "w") as fh:
		for name, seq in asr_seqs.items():
			leaves_str = asr_node_leaves.get(name, "NA")
			fh.write(f">{name}|status={asr_status[name]}|species={leaves_str}\n{seq}\n")



# Per-ORF tree plotting

def _phylo_coords(tree):
	"""Compute x (depth from root) and y coordinates for every clade, matching
	exactly the layout Bio.Phylo.draw uses internally so overlaid dots align
	with the rendered tree."""
	# x: same as Phylo.draw's get_x_positions
	depths = tree.depths()
	if not max(depths.values()):
		depths = tree.depths(unit_branch_lengths=True)
	# y: same as Phylo.draw's get_y_positions — note the reversed() and the
	# postorder "mean of immediate children" (NOT mean of descendant leaves).
	terminals = tree.get_terminals()
	maxheight = len(terminals)
	y_coords = {tip: maxheight - i for i, tip in enumerate(reversed(terminals))}
	for clade in tree.get_nonterminals(order="postorder"):
		if clade.clades:
			y_coords[clade] = (y_coords[clade.clades[0]] + y_coords[clade.clades[-1]]) / 2.0
	return depths, y_coords, maxheight


def _plot_tree_with_dots(newick_str, node_colors, output_path, title,
						 legend_items, show_internal_labels=False):
	"""
	Render a Newick tree with colored dots overlaid on nodes.
	node_colors: {name_or_None: hex_color}  (internal nodes can be included by name)
	legend_items: list of (label, color) tuples for the legend
	"""
	if not _PLOTTING_AVAILABLE:
		return
	tree = Phylo.read(StringIO(newick_str), "newick")
	n_leaves = tree.count_terminals()
	fig_height = max(4.0, 0.13 * n_leaves)
	fig, ax = plt.subplots(figsize=(9, fig_height))

	def label_func(clade):
		if clade.name is None:
			return ""
		if clade.is_terminal():
			return clade.name
		if show_internal_labels and clade.name.startswith("Node"):
			return clade.name
		return ""

	Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
			   label_func=label_func)

	depths, y_coords, _ = _phylo_coords(tree)
	for clade in tree.find_clades():
		name = clade.name
		if name is None:
			continue
		# Biopython's Newick parser may mangle special characters (dots, spaces)
		# in taxon names.  Try the raw name first, then a normalised version.
		color = node_colors.get(name)
		if color is None:
			norm = name.replace("_", " ").replace(".", "_").replace(" ", "_")
			color = node_colors.get(norm)
		if color is None:
			# Last resort: match by stripping non-alphanumeric chars
			stripped = "".join(c if c.isalnum() or c == "_" else "_" for c in name)
			for k, v in node_colors.items():
				k_stripped = "".join(c if c.isalnum() or c == "_" else "_" for c in k)
				if k_stripped == stripped:
					color = v
					break
		if color is None:
			continue
		x = depths.get(clade, 0)
		y = y_coords.get(clade, 0)
		ax.scatter([x], [y], c=[color], s=45, zorder=10,
				   edgecolors="black", linewidths=0.4)

	ax.set_title(title, fontsize=11)
	ax.set_xlabel("")
	ax.set_ylabel("")

	# Legend (custom handles, placed at lower right outside)
	handles = [Line2D([0], [0], marker='o', color='w', label=lbl,
					  markerfacecolor=col, markeredgecolor='black',
					  markeredgewidth=0.4, markersize=8)
			   for lbl, col in legend_items]
	ax.legend(handles=handles, loc="lower right", frameon=True,
			  fontsize=9, framealpha=0.9)

	plt.tight_layout()
	plt.savefig(output_path, dpi=120, bbox_inches="tight")
	plt.close(fig)


def write_orf_trees(asr_tree, full_tree_path, orf_id, outdir, ref_sp,
					species_with_orf_local, asr_fasta,
					stop_cutoff, invalid_cutoff, start_mode,
					initiation_offset, min_identity):
	"""Write two SVGs per ORF: {orf_id}.naive and {orf_id}.asr."""
	if not _PLOTTING_AVAILABLE:
		return

	# Parse species_with_orf_local (semicolon-separated)
	orf_species = set()
	if species_with_orf_local and species_with_orf_local != "NA":
		orf_species = set(s for s in species_with_orf_local.split(";") if s)

	# Identify which species were in the local alignment (from the .anc.fas)
	aligned_species = set()
	try:
		pairs = parse_fasta_str(open(asr_fasta).read())
		aligned_species = {n for n, _ in pairs if not n.startswith("Node")}
	except Exception:
		pass

	# Naive plot: full species tree, leaves colored by ORF status
	try:
		with open(full_tree_path) as fh:
			full_newick = fh.read().strip()
		tfull = Tree(open(full_tree_path), parser=1)
		leaf_names = [l.name for l in tfull.leaves()]
		naive_colors = {}
		for sp in leaf_names:
			if sp in orf_species:
				naive_colors[sp] = COLOR_CONSERVED
			elif sp in aligned_species:
				naive_colors[sp] = COLOR_NOTCONSERVED
			else:
				naive_colors[sp] = COLOR_MISSING
		legend_naive = [
			("ORF conserved",      COLOR_CONSERVED),
			("ORF not conserved",  COLOR_NOTCONSERVED),
			("Not aligned",        COLOR_MISSING),
		]
		_plot_tree_with_dots(
			full_newick, naive_colors,
			os.path.join(outdir, f"{orf_id}.naive.svg"),
			title=f"{orf_id} — naive (species tree)",
			legend_items=legend_naive,
			show_internal_labels=False)
	except Exception as exc:
		LOG.debug("Naive tree plot failed for %s: %s", orf_id, exc)

	# ASR plot: PRANK tree, leaves + Node* colored by orf property
	try:
		with open(asr_tree) as fh:
			asr_newick = fh.read().strip()
		tasr = PhyloTree(open(asr_tree), parser=1, alignment=asr_fasta, alg_format="fasta")
		try:
			ref_seq = tasr[ref_sp].props["sequence"]
		except Exception:
			return
		asr_colors = {}
		for n in tasr.traverse():
			cons = is_conserved(
				ref_seq, n.props.get("sequence", ""),
				stop_cutoff=stop_cutoff, invalid_cutoff=invalid_cutoff,
				start_mode=start_mode, initiation_offset=initiation_offset,
				min_identity=min_identity, query_is_ref=(n.name == ref_sp))
			if not cons.valid:
				asr_colors[n.name] = COLOR_MISSING
			elif cons.conserved:
				asr_colors[n.name] = COLOR_CONSERVED
			else:
				asr_colors[n.name] = COLOR_NOTCONSERVED
		legend_asr = [
			("ORF conserved",      COLOR_CONSERVED),
			("ORF not conserved",  COLOR_NOTCONSERVED),
			("Invalid sequence",   COLOR_MISSING),
		]
		_plot_tree_with_dots(
			asr_newick, asr_colors,
			os.path.join(outdir, f"{orf_id}.asr.svg"),
			title=f"{orf_id} — ASR (PRANK tree)",
			legend_items=legend_asr,
			show_internal_labels=True)
	except Exception as exc:
		LOG.debug("ASR tree plot failed for %s: %s", orf_id, exc)



# Error classification

def _classify_error(exc):
	msg = str(exc)
	if "Could not download" in msg or "URL fetch failed" in msg:
		return "alignment_download_failed"
	if "Timeout" in msg or "timed out" in msg:
		return "timeout"
	if "not in alignment" in msg:
		return "ref_species_not_in_alignment"
	if "Empty alignment" in msg:
		return "empty_alignment"
	if "too short" in msg:
		return "sequence_too_short"
	if "PRANK failed" in msg:
		return "prank_failed"
	if "PRANK timed out" in msg:
		return "prank_timeout"
	if "PRANK did not produce" in msg:
		return "prank_no_output"
	if "CodAlignView error" in msg:
		return "codalignview_error"
	if isinstance(exc, FileNotFoundError):
		return "file_not_found"
	if isinstance(exc, subprocess.TimeoutExpired):
		return "prank_timeout"
	return f"other:{type(exc).__name__}"



# Single-ORF worker

def _process_single_orf(
	orf_id, orf, alnset, ref_species, full_tree_path,
	aln_dir, prank_dir, fasta_dir, trees_dir,
	stop_cutoff, invalid_cutoff, start_mode, initiation_offset,
	min_identity, force, plot_trees,
):
	"""Returns (orf_id, result_dict, error_reason_or_None).
	   result_dict is returned on success and failure."""
	try:
		intervals = gtf_to_intervals(orf)
		aln_fasta_path = os.path.join(aln_dir, f"{orf_id}.fa")
		raw_fasta_path = os.path.join(aln_dir, f"{orf_id}_raw.fa")

		if not os.path.exists(raw_fasta_path) or force:
			fasta_str = download_alignment(intervals, orf["strand"], alnset)
			with open(raw_fasta_path, "w") as fh:
				fh.write(fasta_str)
		else:
			with open(raw_fasta_path, "r") as fh:
				fasta_str = fh.read()

		if not os.path.exists(aln_fasta_path) or force:
			processed_fasta, enough_species = prepare_codon_alignment(
				fasta_str, ref_species=ref_species,
				min_species=2, invalid_cutoff=invalid_cutoff)
			with open(aln_fasta_path, "w") as fh:
				fh.write(processed_fasta)
		else:
			pairs = parse_fasta_str(open(aln_fasta_path).read())
			enough_species = len(pairs) >= 2

		if not enough_species:
			bls = _too_few_species_result(ref_species)
			bls["orf_id"] = orf_id
			bls["chrom"] = orf["chrom"]
			bls["strand"] = orf["strand"]
			bls["n_exons"] = len(orf["exons"])
			return (orf_id, bls, None)

		prank_prefix = os.path.join(prank_dir, orf_id)
		asr_fas = prank_prefix + ".anc.fas"
		asr_dnd = prank_prefix + ".anc.dnd"
		if not (os.path.exists(asr_fas) and os.path.exists(asr_dnd)) or force:
			asr_fas, asr_dnd = run_prank(aln_fasta_path, full_tree_path, prank_prefix)

		bls = compute_bls(
			asr_fasta=asr_fas, asr_tree=asr_dnd, full_tree_path=full_tree_path,
			ref_sp=ref_species, stop_cutoff=stop_cutoff,
			invalid_cutoff=invalid_cutoff, start_mode=start_mode,
			initiation_offset=initiation_offset, min_identity=min_identity)
		bls["orf_id"] = orf_id
		bls["chrom"] = orf["chrom"]
		bls["strand"] = orf["strand"]
		bls["n_exons"] = len(orf["exons"])

		try:
			write_orf_alignments(
				asr_fas, asr_dnd, orf_id, fasta_dir, ref_sp=ref_species,
				stop_cutoff=stop_cutoff, invalid_cutoff=invalid_cutoff,
				start_mode=start_mode, initiation_offset=initiation_offset,
				min_identity=min_identity)
		except Exception:
			pass

		if plot_trees and trees_dir is not None:
			try:
				write_orf_trees(
					asr_dnd, full_tree_path, orf_id, trees_dir, ref_species,
					bls.get("species_with_orf_local", ""), asr_fas,
					stop_cutoff, invalid_cutoff, start_mode,
					initiation_offset, min_identity)
			except Exception as exc:
				LOG.debug("Tree plotting failed for %s: %s", orf_id, exc)

		return (orf_id, bls, None)

	except Exception as exc:
		reason = _classify_error(exc)
		# Failed ORFs still go to main TSV
		bls = _failed_result(ref_species, reason=reason)
		bls["orf_id"] = orf_id
		bls["chrom"] = orf["chrom"]
		bls["strand"] = orf["strand"]
		bls["n_exons"] = len(orf["exons"])
		return (orf_id, bls, reason)



# Main pipeline

COL_ORDER = [
	"orf_id", "chrom", "strand", "n_exons",
	"origin_manner", "origin_age_local", "origin_age_global",
	"n_species_with_orf_local", "species_with_orf_local",
	"most_distant_orf_species_local",
	"identity_conserved", "identity_notconserved", "stop_conserved",
	"repeat_overlaps", "blast_matches",
	"origin_species",
	"bl_all", "bl_sub_origin", "bl_orf_origin",
	"bls_global_origin", "bls_local_origin",
	"origin_age_global_fixation",
	"naive_age_global", "bl_sub_naive", "bl_orf_naive",
	"bls_global_naive", "bls_local_naive",
	"naive_age_global_fixation",
]

def run_pipeline(args):
	outdir = Path(args.outdir)
	outdir.mkdir(parents=True, exist_ok=True)
	aln_dir = str(outdir / "alignments")
	prank_dir = str(outdir / "prank")
	fasta_dir = str(outdir / "fasta")
	trees_dir = str(outdir / "trees") if args.plot_trees else None
	os.makedirs(aln_dir, exist_ok=True)
	os.makedirs(prank_dir, exist_ok=True)
	os.makedirs(fasta_dir, exist_ok=True)
	if trees_dir is not None:
		os.makedirs(trees_dir, exist_ok=True)
		if not _PLOTTING_AVAILABLE:
			LOG.warning("--plot-trees requested but matplotlib/Bio.Phylo not available; skipping plots.")

	orfs = parse_gtf(args.gtf)
	if not orfs:
		LOG.error("No ORFs parsed from GTF.")
		return

	# Load repeat annotations if provided
	repeats = None
	if args.repeats:
		repeats = parse_repeatmasker(args.repeats)

	if args.tree:
		full_tree_path = args.tree
		LOG.info("Using provided species tree: %s", full_tree_path)
	else:
		full_tree_path = str(outdir / "species_tree.nh")
		if os.path.exists(full_tree_path) and not args.force:
			LOG.info("Species tree already exists: %s", full_tree_path)
		else:
			LOG.info("Downloading species tree for alnset '%s' ...", args.alnset)
			tree_str = download_alnset_tree(args.alnset)
			with open(full_tree_path, "w") as fh:
				fh.write(tree_str + "\n")
			LOG.info("Species tree saved to %s", full_tree_path)

	try:
		_t = Tree(open(full_tree_path), parser=1)
		LOG.info("Species tree has %d leaves", len(list(_t.leaves())))
	except Exception as exc:
		LOG.error("Cannot parse species tree %s: %s", full_tree_path, exc)
		return

	n_total = len(orfs)
	threads = max(1, args.threads)
	LOG.info("Processing %d ORFs with %d thread(s)...", n_total, threads)

	results = []
	failed = []
	n_ok = n_fail = 0

	if threads == 1:
		for idx, (orf_id, orf) in enumerate(orfs.items(), 1):
			LOG.info("[%d/%d] Processing %s", idx, n_total, orf_id)
			oid, res, err = _process_single_orf(
				orf_id, orf, args.alnset, args.ref_species,
				full_tree_path, aln_dir, prank_dir, fasta_dir, trees_dir,
				args.stop_cutoff, args.invalid_cutoff,
				args.start_mode, args.initiation_offset,
				args.min_identity, args.force, args.plot_trees)
			results.append(res)
			if err:
				LOG.warning("  SKIP %s: %s", oid, err)
				failed.append((oid, err))
				n_fail += 1
			else:
				n_ok += 1
	else:
		futures = {}
		with ProcessPoolExecutor(max_workers=threads) as pool:
			for orf_id, orf in orfs.items():
				fut = pool.submit(
					_process_single_orf,
					orf_id, orf, args.alnset, args.ref_species,
					full_tree_path, aln_dir, prank_dir, fasta_dir, trees_dir,
					args.stop_cutoff, args.invalid_cutoff,
					args.start_mode, args.initiation_offset,
					args.min_identity, args.force, args.plot_trees)
				futures[fut] = orf_id
			done_count = 0
			for fut in as_completed(futures):
				done_count += 1
				oid = futures[fut]
				try:
					_, res, err = fut.result()
				except Exception as exc:
					err = _classify_error(exc)
					res = _failed_result(args.ref_species, reason=err)
					res["orf_id"] = oid
					res["chrom"] = orfs[oid]["chrom"]
					res["strand"] = orfs[oid]["strand"]
					res["n_exons"] = len(orfs[oid]["exons"])
				results.append(res)
				if err:
					LOG.warning("[%d/%d] SKIP %s: %s", done_count, n_total, oid, err)
					failed.append((oid, err))
					n_fail += 1
				else:
					LOG.info("[%d/%d] Done %s", done_count, n_total, oid)
					n_ok += 1

	# Write ALL results (successful + failed) to main TSV
	# Annotate with repeat overlaps and BLAST
	if results:
		# Repeat overlaps
		for r in results:
			if repeats is not None:
				oid = r["orf_id"]
				if oid in orfs:
					r["repeat_overlaps"] = find_repeat_overlaps(
						orfs[oid]["chrom"], orfs[oid]["exons"], repeats)
				else:
					r["repeat_overlaps"] = "NA"
			else:
				r["repeat_overlaps"] = "NA"

		# BLAST search
		if args.blast:
			LOG.info("Running BLASTP for %d ORFs...", len(results))
			blast_dir = str(outdir / "blast")
			os.makedirs(blast_dir, exist_ok=True)
			blast_db = setup_blast_db(args.blast, str(outdir))
			blast_evalue = args.blast_evalue

			for i, r in enumerate(results):
				oid = r["orf_id"]
				# Get the reference protein sequence from the prot.fa file
				prot_fa = os.path.join(fasta_dir, f"{oid}.prot.fa")
				ref_prot = None
				if os.path.exists(prot_fa):
					pairs = parse_fasta_str(open(prot_fa).read())
					for name, seq in pairs:
						# Find the reference species entry
						if name.startswith(args.ref_species + "|") or name == args.ref_species:
							# Remove stop codon if present
							ref_prot = seq.rstrip("*")
							break
				if ref_prot and len(ref_prot) >= 5:
					r["blast_matches"] = run_blastp(ref_prot, blast_db, blast_evalue,
												   oid, blast_dir)
				else:
					r["blast_matches"] = "none"
				if (i + 1) % 500 == 0:
					LOG.info("  BLAST: %d/%d done", i + 1, len(results))
		else:
			for r in results:
				r["blast_matches"] = "NA"

		results.sort(key=lambda r: r["orf_id"])
		out_tsv = str(outdir / "bls_results.tsv")
		with open(out_tsv, "w") as fh:
			fh.write("\t".join(COL_ORDER) + "\n")
			for r in results:
				fh.write("\t".join(str(r.get(c, "")) for c in COL_ORDER) + "\n")
		LOG.info("Results written to %s (%d rows)", out_tsv, len(results))

	# Also write failed ORFs separately with reasons
	if failed:
		failed.sort(key=lambda x: x[0])
		fail_tsv = str(outdir / "bls_results_noreconstructed.tsv")
		with open(fail_tsv, "w") as fh:
			fh.write("orf_id\treason\n")
			for oid, reason in failed:
				fh.write(f"{oid}\t{reason}\n")
		LOG.info("Failed ORFs written to %s (%d ORFs)", fail_tsv, len(failed))

	LOG.info("Done: %d succeeded, %d failed out of %d ORFs.", n_ok, n_fail, n_total)



# Main

def main():
	parser = argparse.ArgumentParser(
		description="ORF conservation & BLS pipeline: GTF -> CodAlignView -> PRANK -> BLS")
	parser.add_argument("--gtf", required=True,
		help="GTF file with ORFs (orf_id in attributes).")
	parser.add_argument("--alnset", required=True,
		help="CodAlignView alignment set.")
	parser.add_argument("--ref-species", required=True,
		help="Reference species name (e.g. Human).")
	parser.add_argument("--outdir", required=True,
		help="Output directory.")
	parser.add_argument("--tree", default=None,
		help="Species tree (Newick). Auto-downloaded if omitted.")
	parser.add_argument("--stop-cutoff", type=float, default=0.7,
		help="Min fraction of ORF that must be stop-free (default: 0.7).")
	parser.add_argument("--invalid-cutoff", type=float, default=0.5,
		help="Max fraction of N in ungapped sequence (default: 0.5).")
	parser.add_argument("--start-mode", choices=["same", "nearcognate", "atg", "none"],
		default="same",
		help="Start codon criterion: 'same'=match ref or ATG (default), "
			 "'nearcognate'=ATG/CTG/GTG/TTG, 'atg'=ATG only, 'none'=ignored.")
	parser.add_argument("--initiation-offset", type=int, default=3,
		help="Codons from 5' end to search for start codon (default: 3).")
	parser.add_argument("--min-identity", type=float, default=0.0,
		help="Min amino acid identity vs ref to count as conserved (default: 0.0 = off).")
	parser.add_argument("--repeats", default=None,
		help="RepeatMasker .out file. If provided, adds repeat_overlaps column.")
	parser.add_argument("--blast", default=None, metavar="FASTA",
		help="Protein FASTA database for BLASTP. Each ORF's ref protein is searched against it.")
	parser.add_argument("--blast-evalue", type=float, default=1e-4,
		help="E-value threshold for BLAST hits (default: 1e-4).")
	parser.add_argument("--plot-trees", action="store_true",
		help="Generate per-ORF SVG trees in outdir/trees/ (naive + ASR, dots colored by conservation).")
	parser.add_argument("--force", action="store_true",
		help="Re-download and re-run everything.")
	parser.add_argument("--threads", type=int, default=1,
		help="Parallel ORFs (default: 1).")
	parser.add_argument("-v", "--verbose", action="store_true",
		help="Verbose logging.")
	args = parser.parse_args()
	logging.basicConfig(
		level=logging.DEBUG if args.verbose else logging.INFO,
		format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
	run_pipeline(args)

if __name__ == "__main__":
	main()
