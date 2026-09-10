#!/usr/bin/env python3
"""Large-scale backtest: evaluate pipeline back-mutation prediction accuracy
on mouse→humanized antibody pairs across the 9 germline selection strategies.

For each pair:
  - Run 9 germline selection strategies (VH + VL independently)
  - Compare pipeline T1+T2 recommendations vs actual humanized back-mutations
  - Calculate precision/recall/F1 for each strategy
  - Generate summary report with strategy rankings

Data requirement:
  data/benchmarks/humanized_pairs.csv with columns
  name, mouse_vh, mouse_vl, humanized_vh, humanized_vl
  (this benchmark pair table is not bundled with the repository; the script
  exits with a clear message when it is absent).

This is the source of the strategy F1 table in docs/scoring.md section 0.
Run:  python3 tests/backtest_large.py
"""
import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from humanize.backmut import analyze_backmutations
from humanize.germline import GermlineDB, GermlineChoice, load_germline_db
from humanize.graft import graft_chain
from humanize.numbering import number_heavy, number_light
from humanize.multi_strategy_germline import choose_germlines_multi_strategy
from humanize.sequences import InputChain, classify_sequences, detect_format


STRATEGIES = [
    "fr_best", "cdr_best", "composite", "cvi_best",
    "min_backmutations", "current", "adimab_frequency",
    "pioneer_frequency", "composite_3axis",
]


def load_pairs(csv_path):
    """Load mouse→humanized pairs from CSV."""
    pairs = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pairs.append(row)
    return pairs


def analyze_one_chain(parent_seq, actual_seq, chain_type, strategy_name, v_gene, db):
    """Analyze back-mutation prediction accuracy for one chain."""
    try:
        if chain_type == "VH":
            parent = number_heavy(parent_seq)
            actual = number_heavy(actual_seq)
        else:
            parent = number_light(parent_seq)
            actual = number_light(actual_seq)
    except Exception as e:
        return None

    if v_gene is None or v_gene.numbered is None:
        return None

    # Pipeline recommendation on parent
    try:
        bm = analyze_backmutations(parent, v_gene, is_vhh=False)
        t12 = set(bm.revert_positions(("T1", "T2")))
    except Exception:
        return None

    # Actual humanized back-mutations
    dmap = parent.posmap()
    amap = actual.posmap()
    gmap = v_gene.numbered.posmap()

    actual_bm = set()
    tp = fp = fn = tn = 0
    details = []

    for pos in sorted(set(dmap) & set(gmap),
                      key=lambda p: (p[0], int("".join(c for c in p if c.isdigit())))):
        num = int("".join(c for c in pos if c.isdigit()))
        if num >= (103 if chain_type == "VH" else 98):
            continue
        if parent.region_of(pos) not in ("FR1", "FR2", "FR3"):
            continue
        if chain_type == "VH" and num in (93, 94):
            continue
        if dmap[pos] == gmap[pos] or pos not in amap:
            continue

        actual_aa = amap[pos]
        if actual_aa == dmap[pos]:
            is_bm = True
            kind = "donor"
        elif actual_aa == gmap[pos]:
            is_bm = False
            kind = "germline"
        else:
            is_bm = False
            kind = "compromise"

        rec = pos in t12

        if is_bm:
            actual_bm.add(pos)
        if is_bm and rec:
            tp += 1
        elif is_bm and not rec:
            fn += 1
        elif not is_bm and rec:
            fp += 1
        else:
            tn += 1
        details.append((pos, dmap[pos], gmap[pos], actual_aa, kind, rec))

    precision = round(tp / (tp + fp), 3) if tp + fp else None
    recall = round(tp / (tp + fn), 3) if tp + fn else None
    f1 = round(2 * precision * recall / (precision + recall), 3) if precision and recall else None

    return {
        "actual_backmutations": sorted(actual_bm),
        "our_T1T2": sorted(t12),
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "precision": precision, "recall": recall, "f1": f1,
        "details": details,
    }


def backtest_pair(pair, db):
    """Backtest one mouse→humanized pair across all 9 strategies."""
    mouse_vh = pair["mouse_vh"].replace(" ", "")
    mouse_vl = pair["mouse_vl"].replace(" ", "")
    humanized_vh = pair["humanized_vh"].replace(" ", "")
    humanized_vl = pair["humanized_vl"].replace(" ", "")

    results = {s: {"VH": [], "VL": []} for s in STRATEGIES}

    for chain_type, mouse_seq, actual_seq in [
        ("VH", mouse_vh, humanized_vh),
        ("VL", mouse_vl, humanized_vl),
    ]:
        if not mouse_seq or not actual_seq:
            continue
        try:
            if chain_type == "VH":
                parent = number_heavy(mouse_seq)
            else:
                parent = number_light(mouse_seq)
        except Exception:
            continue

        try:
            multi_result = choose_germlines_multi_strategy(parent, db, n_alternatives=1)
        except Exception:
            continue

        for strategy in STRATEGIES:
            candidate = multi_result.get_best(strategy)
            if candidate is None:
                results[strategy][chain_type] = None
                continue
            v_gene = candidate.gene
            res = analyze_one_chain(mouse_seq, actual_seq, chain_type, strategy, v_gene, db)
            results[strategy][chain_type] = res

    return results


def compute_strategy_metrics(strategy_results):
    """Aggregate precision/recall across all pairs for a strategy."""
    all_tp = all_fp = all_fn = all_tn = 0
    chain_metrics = {"VH": {"tp": 0, "fp": 0, "fn": 0, "tn": 0},
                     "VL": {"tp": 0, "fp": 0, "fn": 0, "tn": 0}}
    n_pairs_with_data = 0

    for pair_results in strategy_results:
        for chain_type in ("VH", "VL"):
            res = pair_results[chain_type]
            if res is None:
                continue
            all_tp += res["tp"]
            all_fp += res["fp"]
            all_fn += res["fn"]
            all_tn += res["tn"]
            chain_metrics[chain_type]["tp"] += res["tp"]
            chain_metrics[chain_type]["fp"] += res["fp"]
            chain_metrics[chain_type]["fn"] += res["fn"]
            chain_metrics[chain_type]["tn"] += res["tn"]
        n_pairs_with_data += 1

    precision = round(all_tp / (all_tp + all_fp), 3) if all_tp + all_fp else None
    recall = round(all_tp / (all_tp + all_fn), 3) if all_tp + all_fn else None
    f1 = round(2 * precision * recall / (precision + recall), 3) if precision and recall else None

    vh_precision = round(chain_metrics["VH"]["tp"] / (chain_metrics["VH"]["tp"] + chain_metrics["VH"]["fp"]), 3) if chain_metrics["VH"]["tp"] + chain_metrics["VH"]["fp"] else None
    vh_recall = round(chain_metrics["VH"]["tp"] / (chain_metrics["VH"]["tp"] + chain_metrics["VH"]["fn"]), 3) if chain_metrics["VH"]["tp"] + chain_metrics["VH"]["fn"] else None
    vl_precision = round(chain_metrics["VL"]["tp"] / (chain_metrics["VL"]["tp"] + chain_metrics["VL"]["fp"]), 3) if chain_metrics["VL"]["tp"] + chain_metrics["VL"]["fp"] else None
    vl_recall = round(chain_metrics["VL"]["tp"] / (chain_metrics["VL"]["tp"] + chain_metrics["VL"]["fn"]), 3) if chain_metrics["VL"]["tp"] + chain_metrics["VL"]["fn"] else None

    return {
        "precision": precision, "recall": recall, "f1": f1,
        "tp": all_tp, "fp": all_fp, "fn": all_fn, "tn": all_tn,
        "n_pairs": n_pairs_with_data,
        "vh": {"precision": vh_precision, "recall": vh_recall},
        "vl": {"precision": vl_precision, "recall": vl_recall},
    }


def main():
    csv_path = os.path.join(ROOT, "data", "benchmarks", "humanized_pairs.csv")
    if not os.path.exists(csv_path):
        print(
            "SKIP: benchmark pair table not found:\n"
            f"  {csv_path}\n"
            "This large-scale backtest needs data/benchmarks/humanized_pairs.csv\n"
            "(name, mouse_vh, mouse_vl, humanized_vh, humanized_vl). It is not\n"
            "bundled with the repository. The portable regressions are\n"
            "tests/backtest.py and tests/backtest_scale.py.")
        return 1
    pairs = load_pairs(csv_path)
    print(f"Loaded {len(pairs)} mouse→humanized pairs")

    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    print(f"Loaded germline DB")

    all_pair_results = []
    for i, pair in enumerate(pairs):
        name = pair["name"]
        print(f"  [{i+1}/{len(pairs)}] {name}...")
        try:
            results = backtest_pair(pair, db)
            all_pair_results.append(results)
        except Exception as e:
            print(f"    ERROR: {e}")
            all_pair_results.append({s: {"VH": None, "VL": None} for s in STRATEGIES})

    print(f"\n{'='*80}")
    print("STRATEGY COMPARISON REPORT")
    print(f"{'='*80}\n")

    strategy_metrics = {}
    for strategy in STRATEGIES:
        strategy_results = [pr[strategy] for pr in all_pair_results]
        metrics = compute_strategy_metrics(strategy_results)
        strategy_metrics[strategy] = metrics

    print(f"{'Strategy':<22} {'Precision':>10} {'Recall':>10} {'F1':>10} {'TP':>5} {'FP':>5} {'FN':>5}")
    print("-" * 82)
    for strategy in STRATEGIES:
        m = strategy_metrics[strategy]
        p = f"{m['precision']:.3f}" if m['precision'] is not None else "N/A"
        r = f"{m['recall']:.3f}" if m['recall'] is not None else "N/A"
        f = f"{m['f1']:.3f}" if m['f1'] is not None else "N/A"
        print(f"{strategy:<22} {p:>10} {r:>10} {f:>10} {m['tp']:>5} {m['fp']:>5} {m['fn']:>5}")

    print(f"\n\n{'='*80}")
    print("CHAIN-SPECIFIC BREAKDOWN")
    print(f"{'='*80}\n")

    print(f"{'Strategy':<22} {'VH Prec':>10} {'VH Rec':>10} {'VL Prec':>10} {'VL Rec':>10}")
    print("-" * 72)
    for strategy in STRATEGIES:
        m = strategy_metrics[strategy]
        vhp = f"{m['vh']['precision']:.3f}" if m['vh']['precision'] is not None else "N/A"
        vhr = f"{m['vh']['recall']:.3f}" if m['vh']['recall'] is not None else "N/A"
        vlp = f"{m['vl']['precision']:.3f}" if m['vl']['precision'] is not None else "N/A"
        vlr = f"{m['vl']['recall']:.3f}" if m['vl']['recall'] is not None else "N/A"
        print(f"{strategy:<22} {vhp:>10} {vhr:>10} {vlp:>10} {vlr:>10}")

    ranked = sorted(strategy_metrics.items(), key=lambda x: x[1]['f1'] or 0, reverse=True)
    print(f"\n\n{'='*80}")
    print("BEST STRATEGY (by F1)")
    print(f"{'='*80}\n")
    for rank, (strategy, m) in enumerate(ranked, 1):
        print(f"  #{rank} {strategy}: F1={m['f1']:.3f}, P={m['precision']:.3f}, R={m['recall']:.3f}")

    print(f"\n\n{'='*80}")
    print("PER-PAIR DETAILS")
    print(f"{'='*80}\n")

    for i, (pair, pair_results) in enumerate(zip(pairs, all_pair_results)):
        name = pair["name"]
        print(f"[{i+1}] {name}")
        for strategy in STRATEGIES:
            sr = pair_results[strategy]
            vh_res = sr["VH"]
            vl_res = sr["VL"]
            vh_info = f"VH: T1+T2={vh_res['our_T1T2']}, actual={vh_res['actual_backmutations']}, P={vh_res['precision']}, R={vh_res['recall']}" if vh_res else "VH: N/A"
            vl_info = f"VL: T1+T2={vl_res['our_T1T2']}, actual={vl_res['actual_backmutations']}, P={vl_res['precision']}, R={vl_res['recall']}" if vl_res else "VL: N/A"
            print(f"  {strategy:<22} {vh_info}")
            print(f"  {'':<22} {vl_info}")
        print()

    out_path = os.path.join(ROOT, "outputs", "backtest_large_results.txt")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write("LARGE-SCALE BACKTEST RESULTS\n")
        f.write(f"Pairs: {len(pairs)}, Strategies: {len(STRATEGIES)}\n\n")
        f.write(f"{'Strategy':<22} {'Precision':>10} {'Recall':>10} {'F1':>10} {'TP':>5} {'FP':>5} {'FN':>5}\n")
        f.write("-" * 82 + "\n")
        for strategy in STRATEGIES:
            m = strategy_metrics[strategy]
            p = f"{m['precision']:.3f}" if m['precision'] is not None else "N/A"
            r = f"{m['recall']:.3f}" if m['recall'] is not None else "N/A"
            f1_str = f"{m['f1']:.3f}" if m['f1'] else "N/A"
            f.write(f"{strategy:<22} {p:>10} {r:>10} {f1_str:>10} {m['tp']:>5} {m['fp']:>5} {m['fn']:>5}\n")
        f.write(f"\n\nBest by F1: {ranked[0][0]} (F1={ranked[0][1]['f1']:.3f})\n")
    print(f"Results saved to {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
