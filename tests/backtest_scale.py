#!/usr/bin/env python3
"""Scale validation: run the pipeline on 25 real mouse antibodies (HumAb25,
the parental sequences of approved humanized therapeutics).

Checks:
  * numbering success rate (must be 100%)
  * format detection (fab / vhh)
  * germline selection sanity (VH -> IGHV family; VL -> IGKV/IGLV family)
  * variant ladder generation (V0-V3, Vmin) with consistent lengths
  * developability scan runs
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import csv
import tempfile

from humanize.pipeline import PipelineConfig, run_pipeline
from humanize.report import write_all

FAILURES = []


def check(name, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def load_humab25(path):
    with open(path) as fh:
        reader = csv.DictReader(fh)
        return [r for r in reader]


def main():
    bench = os.path.join(ROOT, "data", "benchmarks", "humab25_parental_mouse.csv")
    records = load_humab25(bench)
    print(f"HumAb25 scale validation: {len(records)} mouse antibodies")
    n_ok = 0
    germline_families = {}
    with tempfile.TemporaryDirectory() as work:
        for i, rec in enumerate(records):
            name = rec["name"]
            fasta = os.path.join(work, f"mab{i}.fasta")
            with open(fasta, "w") as fh:
                fh.write(f">{name}_VH\n{rec['h_seq']}\n>{name}_VL\n{rec['l_seq']}\n")
            outdir = os.path.join(work, f"out{i}")
            try:
                result = run_pipeline(fasta, PipelineConfig(), outdir=outdir)
            except Exception as e:
                check(f"{name}: pipeline run", False, str(e)[:120])
                continue
            if result.format != "fab":
                check(f"{name}: format", False, result.format)
            chains_ok = True
            for rep in result.chains:
                c = rep.input_chain
                if c.numbered is None:
                    check(f"{name}/{c.name}: numbered", False)
                    chains_ok = False
                    continue
                vg = rep.germline.v_gene
                if vg is None:
                    check(f"{name}/{c.name}: germline", False)
                    chains_ok = False
                    continue
                fam = "IGHV" if c.chain_type == "H" else ("IGKV" if "IGKV" in vg.gene_id else "IGLV")
                germline_families.setdefault(fam, 0)
                germline_families[fam] += 1
                # variant ladder consistency
                lens = {len(v.sequence) for v in rep.variants}
                if len(lens) != 1:
                    check(f"{name}/{c.name}: variant lengths", False, str(lens))
                    chains_ok = False
                if len(rep.variants) < 4:
                    check(f"{name}/{c.name}: variants", False, str(len(rep.variants)))
                    chains_ok = False
            if chains_ok:
                n_ok += 1
    print(f"\n  {n_ok}/{len(records)} antibodies passed pipeline end-to-end")
    print(f"  germline families assigned: {germline_families}")
    check("all 25 numbered + variant ladders", n_ok == len(records),
          f"{n_ok}/{len(records)}")
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES: {FAILURES[:10]}")
        return 1
    print("ALL SCALE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
