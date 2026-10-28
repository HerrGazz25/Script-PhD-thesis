#!/usr/bin/env python3
"""
Confronto per sequenza tra l'output dello script CDR3 (csv *_cdr3_aa_counts.csv)
e l'export TSV di Platforma, su un campione target e un campione controllo.

Nei TSV di Platforma la CDR3 e' di 21 aa: 'CAK' + 14 aa variabili + 'FDYW'.
Lo script estrae solo i 14 aa variabili, quindi qui i flanking vengono tagliati
(--trim-left 3, --trim-right 4) e le CDR3 con lunghezza diversa da 21 aa o non
produttive vengono escluse, per rendere confrontabili i due strumenti.

Uso:
    python confronta_script_platforma.py \
        --script-target results/target_cdr3_aa_counts.csv \
        --script-control results/controllo_cdr3_aa_counts.csv \
        --plat-target platforma_target.tsv \
        --plat-control platforma_controllo.tsv \
        --out confronto_p3
"""

import argparse

import numpy as np
import pandas as pd


def load_script(path):
    df = pd.read_csv(path)
    return df.groupby("cdr3_aa")["count"].sum()


def load_platforma(path, trim_left, trim_right, expected_len):
    df = pd.read_csv(path, sep="\t",
                     usecols=["readCount", "aaSeqCDR3", "aaLengthCDR3", "isProductiveCDR3"])
    n_all = int(df["readCount"].sum())
    df = df[df["isProductiveCDR3"].astype(str).str.lower() == "true"]
    if expected_len:
        df = df[df["aaLengthCDR3"] == expected_len]
    end = -trim_right if trim_right > 0 else None
    df = df.assign(core=df["aaSeqCDR3"].str[trim_left:end])
    counts = df.groupby("core")["readCount"].sum()
    print(f"    [INFO] {path}: {n_all} reads totali, {int(counts.sum())} dopo i filtri, "
          f"{len(counts)} sequenze aa uniche")
    return counts


def classify(t, c, min_t, min_expected, min_log2fc):
    m = pd.concat([t.rename("count_t"), c.rename("count_c")], axis=1).fillna(0)
    nt, nc = m["count_t"].sum(), m["count_c"].sum()
    m["log2fc"] = np.log2(((m["count_t"] + 0.5) / nt) / ((m["count_c"] + 0.5) / nc))
    m["atteso_c"] = m["count_t"] / nt * nc
    conds = [
        m["count_t"] < min_t,
        (m["count_c"] == 0) & (m["atteso_c"] >= min_expected),
        (m["count_c"] == 0),
        m["log2fc"] >= min_log2fc,
    ]
    labels = ["bassa_copertura", "specifica", "assenza_non_informativa", "arricchita"]
    m["classe"] = np.select(conds, labels, default="altro")
    return m


def compare_counts(name, s, p, min_count):
    m = pd.concat([s.rename("script"), p.rename("platforma")], axis=1).fillna(0)
    both = m[(m["script"] > 0) & (m["platforma"] > 0)]
    solo_s = int(((m["script"] > 0) & (m["platforma"] == 0)).sum())
    solo_p = int(((m["script"] == 0) & (m["platforma"] > 0)).sum())
    top = m[(m["script"] >= min_count) | (m["platforma"] >= min_count)]
    print(f"  [{name}] reads: script={int(m['script'].sum())} platforma={int(m['platforma'].sum())}")
    print(f"  [{name}] sequenze in comune={len(both)}, solo script={solo_s}, solo platforma={solo_p}")
    print(f"  [{name}] Spearman sui conteggi (sequenze con >= {min_count} reads in uno dei due): "
          f"{top['script'].corr(top['platforma'], method='spearman'):.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--script-target", required=True)
    ap.add_argument("--script-control", required=True)
    ap.add_argument("--plat-target", required=True)
    ap.add_argument("--plat-control", required=True)
    ap.add_argument("--out", required=True, help="prefisso dei file di output")
    ap.add_argument("--trim-left", type=int, default=3)
    ap.add_argument("--trim-right", type=int, default=4)
    ap.add_argument("--expected-len", type=int, default=21,
                    help="lunghezza attesa della CDR3 Platforma in aa (0 = non filtrare)")
    ap.add_argument("--min-count", type=int, default=10,
                    help="soglia reads nel target e per il confronto dei conteggi")
    ap.add_argument("--min-expected", type=float, default=5.0,
                    help="reads attese nel controllo per considerare informativa un'assenza")
    ap.add_argument("--min-log2fc", type=float, default=2.0)
    args = ap.parse_args()

    print("Caricamento file...")
    s_t = load_script(args.script_target)
    s_c = load_script(args.script_control)
    p_t = load_platforma(args.plat_target, args.trim_left, args.trim_right, args.expected_len)
    p_c = load_platforma(args.plat_control, args.trim_left, args.trim_right, args.expected_len)

    print("\nConfronto dei conteggi per sequenza:")
    compare_counts("target", s_t, p_t, args.min_count)
    compare_counts("controllo", s_c, p_c, args.min_count)

    cls_s = classify(s_t, s_c, args.min_count, args.min_expected, args.min_log2fc)
    cls_p = classify(p_t, p_c, args.min_count, args.min_expected, args.min_log2fc)
    out = cls_s.add_suffix("_script").join(cls_p.add_suffix("_platforma"), how="outer")
    out["classe_script"] = out["classe_script"].fillna("assente")
    out["classe_platforma"] = out["classe_platforma"].fillna("assente")
    out.index.name = "cdr3_aa"
    out.sort_values("log2fc_platforma", ascending=False).to_csv(f"{args.out}_confronto.csv")

    print("\nClassificazione (righe = script, colonne = Platforma):")
    print(pd.crosstab(out["classe_script"], out["classe_platforma"]))
    print(f"\nTabella completa: {args.out}_confronto.csv")


if __name__ == "__main__":
    main()
