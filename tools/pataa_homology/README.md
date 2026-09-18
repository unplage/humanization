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
| After CDR-quality filter + dedup | **122,610** |
| — VH / VL | 75,448 / 47,162 |
| Top V germlines | IGHV3-23\*01 (16,625), IGKV1-39\*01 (6,541), IGHV1-46\*01 (5,986) |

## Pipeline

`preprocess_pataa.py` — all stages resumable (delete a stage output to re-run):

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
```

Measured on 112 cores (wall time):

| Stage | Time | Notes |
|---|---|---|
| dump | 79 s | 3.98 M seqs, 1.2 GB FASTA |
| split | 11 s | 200 chunks |
| detect | 357 s | 64 workers, `-outfmt 7`, bitscore >= 40 |
| extract | 5 s | 291,360 sequences |
| number | 345 s | 64 workers, ANARCI |
| build | 22 s | filter + dedup |

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

Required tooling (already present in this repo): igblast 1.22.0
(`data/igblast/`), ANARCI, BioPython. `IGDATA` is set automatically.

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
├── preprocess_pataa.py   # dump/split/detect/extract/number/build
├── query_pataa.py        # CDR homology query + reports
└── database/             # pataa.tar.gz + extracted BLAST files (gitignored)
```
