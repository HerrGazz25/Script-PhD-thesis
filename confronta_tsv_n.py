#!/usr/bin/env python3
"""
Confronta due file .tsv:
  - colonna 1 = numero di reads
  - colonna 3 = sequenza

Output (nella cartella di output):
  - riepilogo.csv            numeri e percentuali
  - venn_sequenze.png        Venn sul numero di sequenze uniche (distinte)
  - proporzioni_reads.png    % di reads in sequenze condivise vs esclusive
  - sequenze_condivise.tsv   sequenze presenti in entrambi, con reads dei due file
  - solo_A.tsv / solo_B.tsv  sequenze esclusive di ciascun file

Uso:
  python confronta_tsv.py file_A.tsv file_B.tsv
  python confronta_tsv.py file_A.tsv file_B.tsv --header --nomeA Libreria1 --nomeB Libreria2 --out risultati

Dipendenze:
  pip install pandas matplotlib matplotlib-venn
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
from matplotlib_venn import venn2, venn2_circles


def carica(path, header):
    """Legge solo colonna 1 e 3 e somma le reads per sequenza identica."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=0 if header else None,
        usecols=[0, 2],  # prima e terza colonna
    )
    df.columns = ["reads", "seq"]  # usecols restituisce le colonne in ordine
    df = df.dropna(subset=["seq"])
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0)
    # se la stessa sequenza compare in piu' righe, somma le reads
    return df.groupby("seq", sort=False)["reads"].sum()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file_a")
    p.add_argument("file_b")
    p.add_argument("--header", action="store_true", help="i file hanno una riga di intestazione")
    p.add_argument("--nomeA", default="File A")
    p.add_argument("--nomeB", default="File B")
    p.add_argument("--out", default="output_confronto")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    nA, nB = args.nomeA, args.nomeB

    a = carica(args.file_a, args.header)
    b = carica(args.file_b, args.header)

    idx_shared = a.index.intersection(b.index)
    idx_only_a = a.index.difference(b.index)
    idx_only_b = b.index.difference(a.index)

    # ---- conteggi sequenze distinte ----
    n_only_a, n_only_b, n_shared = len(idx_only_a), len(idx_only_b), len(idx_shared)

    # ---- reads ----
    tot_a, tot_b = a.sum(), b.sum()
    reads_a_shared = a.loc[idx_shared].sum()
    reads_a_only = a.loc[idx_only_a].sum()
    reads_b_shared = b.loc[idx_shared].sum()
    reads_b_only = b.loc[idx_only_b].sum()

    def pct(x, tot):
        return 100 * x / tot if tot else 0

    riepilogo = pd.DataFrame(
        [
            [nA, len(a), n_only_a, n_shared, pct(n_only_a, len(a)), pct(n_shared, len(a)),
             tot_a, reads_a_only, reads_a_shared, pct(reads_a_only, tot_a), pct(reads_a_shared, tot_a)],
            [nB, len(b), n_only_b, n_shared, pct(n_only_b, len(b)), pct(n_shared, len(b)),
             tot_b, reads_b_only, reads_b_shared, pct(reads_b_only, tot_b), pct(reads_b_shared, tot_b)],
        ],
        columns=[
            "file", "seq_distinte", "seq_esclusive", "seq_condivise",
            "%seq_esclusive", "%seq_condivise",
            "reads_totali", "reads_in_seq_esclusive", "reads_in_seq_condivise",
            "%reads_in_seq_esclusive", "%reads_in_seq_condivise",
        ],
    )
    riepilogo.to_csv(os.path.join(args.out, "riepilogo.csv"), index=False)
    print(riepilogo.T.to_string(header=False))

    # ---- Venn sulle sequenze distinte ----
    # Venn a 2 insiemi con aree proporzionali (cerchi e sovrapposizione in scala),
    # stile: cerchi semitrasparenti con bordo nero, legenda a destra,
    # numero di sequenze comuni sotto i cerchi.
    def fmt(n):
        return f"{n:,}".replace(",", ".")  # separatore migliaia all'italiana

    COL_A, COL_B = "#5FD8A4", "#FFA566"  # verde a sinistra, arancione a destra
    subsets = (n_only_a, n_only_b, n_shared)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    v = venn2(subsets=subsets, set_labels=("", ""), set_colors=(COL_A, COL_B), alpha=0.55, ax=ax)
    circ_a, circ_b = venn2_circles(subsets=subsets, ax=ax, color="black", linewidth=1.3)

    # geometria dei cerchi (coordinate dati)
    (xa, ya), ra = circ_a.center, circ_a.radius
    (xb, yb), rb = circ_b.center, circ_b.radius
    span_y = max(ya + ra, yb + rb) - min(ya - ra, yb - rb)

    # numero esclusivo di A: dentro il cerchio verde (posizione di default)
    lab_a = v.get_label_by_id("10")
    if lab_a is not None:
        lab_a.set_text(fmt(n_only_a))
        lab_a.set_fontsize(16)

    # numero esclusivo di B: FUORI dal cerchio arancione, a destra, alla stessa altezza
    lab_b = v.get_label_by_id("01")
    if lab_b is not None:
        lab_b.set_text(fmt(n_only_b))
        lab_b.set_fontsize(16)
        lab_b.set_position((xb + rb + 0.03 * span_y, yb))
        lab_b.set_ha("left")
        lab_b.set_va("center")

    # intersezione: nessun numero all'interno
    lab_ab = v.get_label_by_id("11")
    if lab_ab is not None:
        lab_ab.set_text("")

    # numero di sequenze comuni: sotto i cerchi, centrato orizzontalmente sull'intersezione
    if n_shared > 0 and (xa + ra) > (xb - rb):
        x_inter = ((xb - rb) + (xa + ra)) / 2  # centro della "lente" di sovrapposizione
    else:
        x_inter = (xa + xb) / 2  # nessuna sovrapposizione: punto medio tra i cerchi
    y_sotto = min(ya - ra, yb - rb) - 0.05 * span_y
    ax.text(x_inter, y_sotto, f"{fmt(n_shared)}\nsequenze comuni",
            ha="center", va="top", multialignment="center", fontsize=16)

    # legenda a destra
    handles = [
        Patch(facecolor=COL_A, edgecolor="black", alpha=0.55, label=nA),
        Patch(facecolor=COL_B, edgecolor="black", alpha=0.55, label=nB),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
              frameon=False, fontsize=13)

    fig.savefig(os.path.join(args.out, "venn_sequenze.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---- Barre impilate: % reads in sequenze esclusive vs condivise ----
    fig, ax = plt.subplots(figsize=(5, 5))
    nomi = [nA, nB]
    escl = [pct(reads_a_only, tot_a), pct(reads_b_only, tot_b)]
    cond = [pct(reads_a_shared, tot_a), pct(reads_b_shared, tot_b)]
    ax.bar(nomi, escl, label="Sequenze esclusive")
    ax.bar(nomi, cond, bottom=escl, label="Sequenze condivise")
    ax.set_ylabel("% delle reads totali")
    ax.set_ylim(0, 100)
    ax.set_title("Proporzione di reads per tipo di sequenza")
    ax.legend()
    fig.savefig(os.path.join(args.out, "proporzioni_reads.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---- liste di sequenze ----
    shared_df = pd.DataFrame({f"reads_{nA}": a.loc[idx_shared], f"reads_{nB}": b.loc[idx_shared]})
    shared_df.index.name = "seq"
    shared_df.sort_values(f"reads_{nA}", ascending=False).to_csv(
        os.path.join(args.out, "sequenze_condivise.tsv"), sep="\t")
    a.loc[idx_only_a].sort_values(ascending=False).to_csv(
        os.path.join(args.out, "solo_A.tsv"), sep="\t", header=["reads"])
    b.loc[idx_only_b].sort_values(ascending=False).to_csv(
        os.path.join(args.out, "solo_B.tsv"), sep="\t", header=["reads"])

    print(f"\nFatto. Risultati in: {args.out}/")


if __name__ == "__main__":
    main()
