# Germline data sources

- `human_germline_kabat.json` — **primary bundled source.** 373 human V genes
  + 28 J genes (IGHV/IGKV/IGLV/IGHJ/IGKJ/IGLJ), numbered in Kabat space by
  this pipeline's engine. Derived from the human IMGT germline set embedded in
  the `abnumber` package (MIT license, https://github.com/prihoda/AbNumber),
  which is the same IMGT data used by IgBLAST's `human_gl_*.fasta` files.
  Regenerate with `python3 scripts/build_germline_kabat.py`.

- `abnumber_human_imgt.json` — raw extraction from abnumber (IMGT positions),
  kept for traceability and regeneration.

- On the server, `humanize setup-germline` downloads the NCBI IgBLAST FASTA
  files (`human_gl_V.fasta`, `human_gl_L.fasta`, `human_gl_J.fasta`) from
  `ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/LATEST/`; the
  pipeline prefers those when present and falls back to the bundled set.
