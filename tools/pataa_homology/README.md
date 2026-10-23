# PATAA Patent Antibody Homology Tool

Mine the GenBank Patent Division protein database (**PATAA**, ~4.0 M sequences)
for antibody V domains, build a CDR reference database, and analyse candidate
antibody sequences against it — mirroring the `tools/cdr_homology` workflow but
against the **patent** (IP/engineering) landscape instead of clinical antibodies.

## Database identity

| | |
|---|---|
| Source | GenBank Patent Division (`PATAA`) BLAST v5 |
| Archive | `database/pataa.tar.gz` -> `database/pataa.*` |
| Scale | 3,982,722 sequences; 787,058,342 residues |
| Update | 2026-09-14 |
| Taxonomy | `taxonomy4blast.sqlite3` (6,725 taxids) + `taxdb.btd/.bti` |
| Descriptions | `.phr` deflines, e.g. `Sequence 20020 from patent US 8268964` |

## Built database

`data/pataa_cdrs.json` (+ `data/pataa_stats.json`), **not committed** (regenerate):

| Metric | Value |
|---|---|
| Positive hits (igblastp V bitscore >= 40, human germline) | 289,934 |
| ANARCI-numbered records | 271,114 |
| After CDR-quality filter + dedup | **122,540** |
| — VH / VL | 75,378 / 47,162 |
| Top V germlines | IGHV3-23\*01 (16,423), IGKV1-39\*01 (6,487), IGHV1-46\*01 (5,864) |

The quality filter requires CDR1 ≥ 4, CDR2 ≥ 2 and CDR3 4–40 aa. From the
271,114 numbered records it drops 3,891 with too-short CDRs and 140 with
CDR3 > 40 aa (mostly toxin/cytokine fusions whose C-terminal partner ANARCI
numbers as 111/112 insertions), then dedups by CDR triple.

## Pipeline

`preprocess_pataa.py` — stages are resumable and write outputs atomically
(delete a stage output to re-run):

```bash
# full build (~15 min total on a many-core machine)
python3 tools/pataa_homology/preprocess_pataa.py all --workers 64

# or stage by stage
python3 tools/pataa_homology/preprocess_pataa.py dump      # blastdbcmd -> pataa.fasta (~1.3 min)
python3 tools/pataa_homology/preprocess_pataa.py split     # unique ids -> chunks/ + index.tsv
python3 tools/pataa_homology/preprocess_pataa.py detect    # igblastp, human_V, outfmt 7
python3 tools/pataa_homology/preprocess_pataa.py extract   # positive sequences
python3 tools/pataa_homology/preprocess_pataa.py number    # ANARCI IMGT
python3 tools/pataa_homology/preprocess_pataa.py build     # dedup -> data/pataa_cdrs.json

# after changing CDR extraction: re-extract from cdr_records.jsonl, then build
python3 tools/pataa_homology/preprocess_pataa.py repair
```

`build` trusts the `cdrs` already stored in `cdr_records.jsonl`; `number`
skips when that file exists. So a CDR-extraction fix is **not** picked up by
`build` alone. `repair` streams the JSONL, recomputes `cdrs` from the stored
`numbered` ANARCI posmaps in sequence order, atomically replaces the JSONL,
and runs `build` — no ANARCI or source BLAST DB needed (~1–2 min). It is not
part of `all`.

### Fresh machine rebuild

`data/igblast/` and `database/pataa.*` are **gitignored**, so a clean clone
does not contain them:

```bash
# 1. ANARCI + HMMER + BioPython (or: pip install anarci biopython + hmmer 3.x)
bash scripts/install.sh --full

# 2. IgBLAST 1.22.0 (bin/ + internal_data/); either extract under data/igblast/
#    https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/LATEST/
mkdir -p data/igblast
tar xzf ncbi-igblast-1.22.0-x64-linux.tar.gz -C data/igblast
# ...or point at an existing install:
export IGBLAST_HOME=/path/to/ncbi-igblast-1.22.0

# 3. PATAA archive; dump auto-extracts pataa.* + taxdb.*
cp /path/to/pataa.tar.gz tools/pataa_homology/database/

# 4. full rebuild
python3 tools/pataa_homology/preprocess_pataa.py all --workers 64
```

`IGBLAST_HOME` may be replaced by `IGDATA`, `BLASTDBCMD`, `IGBLASTP` and
`GERMLINE_V` individually. A fresh `all` run uses the fixed ANARCI-order CDR
extraction in `number`, so **`repair` is only needed when reusing a
`cdr_records.jsonl` produced by an older version** (or delete it and re-run
`number`). The expected output is ~122,540 entries; `VH_CDR3` is capped at 40 aa.

Measured on 112 cores (wall time):

| Stage | Time | Notes |
|---|---|---|
| dump | 79 s | 3.98 M seqs, 1.2 GB FASTA |
| split | 11 s | 200 chunks |
| detect | 357 s | 64 workers, `-outfmt 7`, bitscore >= 40 |
| extract | 5 s | 291,360 sequences |
| number | 345 s | 64 workers, ANARCI |
| build | 22 s | CDR-length filter + dedup |

### Why this pipeline

- **`blastdbcmd` first** turns the binary BLAST DB into plain FASTA in ~1 min.
- **`igblastp` is the fast detector** (~345 seq/s/process). It *recognizes the
  Ig domain across species* — mouse 4D5 and camelid VHH are both detected
  against the **human** germline library, so human-only detection is used per
  the requested design. `-num_threads` does **not** help (BLAST parallelises
  within a query, not across a tiny germline DB), so detection is chunked and
  run with **process-level parallelism**.
- **ANARCI does numbering/CDR extraction** because igblastp's germline
  alignment does not yield the full CDR3 (germline lacks CDR3).
- **Dedup by CDR triple** merges the heavy patent redundancy.

Required tooling: IgBLAST 1.22.0 (`data/igblast/` or `IGBLAST_HOME`),
ANARCI + HMMER, BioPython. `IGDATA` is set automatically; see
**Fresh machine rebuild** above.

## Query

```bash
# single pair
python3 tools/pataa_homology/query_pataa.py \
    --vh EVQLQQSGPELVKPGASVKMSCKAS... \
    --vl DIQMTQTTSSLSASLGDRVTISCRAS... \
    --top 5 --output outputs/pataa_homology/

# from FASTA (VH/VL paired by header suffix)
python3 tools/pataa_homology/query_pataa.py \
    --input data/pataa_work/amg110_paired.fasta --top 5 \
    --output outputs/pataa_homology_amg110/
```

Outputs: `pataa_homology_report.md`, `pataa_homology_report.docx`,
`pataa_homology_results.json`, `detail/*.csv`. Each CDR (IMGT) is aligned
(BLOSUM62) to database CDRs; hits are reported with accession, **patent
number(s)**, V-germline and identity. The Word report mirrors
`tools/cdr_homology/query.py` (summary table, per-query CDR tables, top-N per
CDR, appendix). Pass `--no-docx` to skip it.

This module is **standalone** — it is intentionally not wired into the
`scripts/humanize` pipeline.

Implementation note: a k-mer (3-mer) + length prefilter selects the top ~400
candidates per query CDR before full BLOSUM62 alignment, so a whole database
scan is ~3 s per query instead of minutes.

## Interpreting results / caveats

- A **100% CDR identity hit against a patent** is an IP/FTO signal, not proof of
  infringement; review the claim scope before drawing conclusions.
- `species` is blank for most entries (patent deflines often omit it).
- Detection threshold is bitscore >= 40 (bimodal: noise median ~24, real
  antibodies p90 ~179). Raising it trades recall for precision.
- Human-only germline detection was used by design; it detects xenogeneic V
  domains but germline *assignment* for non-human donors may be approximate.

## Files

```
tools/pataa_homology/
├── README.md
├── preprocess_pataa.py   # dump/split/detect/extract/number/repair/build
├── query_pataa.py        # CDR homology query + reports
└── database/             # pataa.tar.gz + extracted BLAST files (gitignored)
```
