#!/usr/bin/env python3
"""
Pipeline completa per library phage display (scFv fusa a pIII), Illumina
NextSeq 500/550 (chimica a 2 colori), indici QIAseq UDI Y-Adapter, amplicone
pre-indici di 196 bp sequenziato paired-end 150x2.

Due modalita' d'uso:

  A) BATCH multi-campione (consigliata per piu' campioni multi-lane):
     legge una cartella con i fastq grezzi Illumina standard
       <codice>_S<n>_L00<lane>_R1_001.fastq.gz  (e relativi R2)
     per ciascun campione: concatena le lane, fa il merge R1+R2 con fastp
     (sfruttando l'overlap sulla CDR3), poi esegue l'estrazione/conteggio
     delle CDR3 sul fastq merged. Un file di riepilogo aggregato viene
     scritto in fondo.

       python phage_cdr3_full_pipeline.py batch \
           --raw-dir raw_fastq \
           --work-dir work \
           --results-dir results \
           --reference pIII_scFv_reference.fasta

  B) SINGOLO file gia' pronto (un fastq unico, gia' mergiato/preparato):

       python phage_cdr3_full_pipeline.py single \
           --fastq sample_merged.fastq.gz \
           --reference pIII_scFv_reference.fasta \
           --out-prefix results/sample1

Dipendenze:
    pip install biopython
    pip install regex          # opzionale ma consigliato per velocita'
    fastp installato e nel PATH # solo per la modalita' batch (merge R1+R2)

Assunzioni sulla CDR3 (confermate dal disegno della library):
 - lunghezza COSTANTE, unico segmento variabile dell'inserto
 - posizioni di default nel riferimento: nt 361-402 (42 nt = 14 aa)
 - le regioni immediatamente a monte/valle sono costanti e vengono usate
   come "anchor" per individuare la CDR3 in ciascuna read, tollerando
   qualche mismatch da errore di sequenziamento

Note NextSeq (chimica a 2 colori):
 - i cicli senza segnale vengono chiamati 'G' con qualita' Phred spesso
   ingannevolmente alta; lo script rimuove le code di poli-G PRIMA del
   quality trimming a finestra scorrevole (--poly-g-min-run, default 10)
 - con amplicone di 196 bp e reads da 150 bp, R1 copre nt 1-150 e R2 copre
   (in coordinate dell'inserto) nt 47-196: la CDR3 (361-402) e' quindi letta
   due volte indipendentemente nell'overlap, da cui il merge con fastp
"""

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

try:
    import regex
    HAVE_REGEX = True
except ImportError:
    HAVE_REGEX = False


# ==========================================================================
# Utility di I/O
# ==========================================================================

def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_reference(fasta_path):
    with open_maybe_gzip(fasta_path) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        sys.exit(f"Nessuna sequenza trovata in {fasta_path}")
    if len(records) > 1:
        print(f"[WARN] Trovate {len(records)} sequenze nel file di riferimento, "
              f"uso la prima ({records[0].id}).", file=sys.stderr)
    return str(records[0].seq).upper()


def count_fastq_reads(path):
    """Conta le read in un fastq (anche .gz) contando le righe / 4."""
    with open_maybe_gzip(path) as handle:
        n_lines = sum(1 for _ in handle)
    return n_lines // 4


# ==========================================================================
# Trimming code di poli-G (artefatto chimica a 2 colori: NextSeq/NovaSeq/iSeq)
# ==========================================================================

def build_poly_g_trimmer(min_run=10):
    if min_run <= 0:
        return lambda record: record
    pattern = re.compile(f"G{{{min_run},}}$")

    def trim(record):
        seq_str = str(record.seq)
        m = pattern.search(seq_str)
        if m:
            return record[:m.start()]
        return record
    return trim


# ==========================================================================
# Quality trimming (sliding window, stile Trimmomatic)
# ==========================================================================

def sliding_window_trim(record, window=4, min_qual=20, min_len=50):
    quals = record.letter_annotations["phred_quality"]
    n = len(quals)
    if n == 0:
        return None

    cut = n
    for i in range(0, n - window + 1):
        avg = sum(quals[i:i + window]) / window
        if avg < min_qual:
            cut = i
            break

    trimmed = record[:cut]
    if len(trimmed) < min_len:
        return None
    return trimmed


# ==========================================================================
# Ricerca approssimata dell'anchor (con mismatch tollerate)
# ==========================================================================

def build_anchor_matcher(anchor, max_mismatches):
    if max_mismatches == 0:
        def find_exact(read_seq):
            idx = read_seq.find(anchor)
            return idx if idx != -1 else None
        return find_exact

    if HAVE_REGEX:
        # s<=N: solo sostituzioni (distanza di Hamming), coerente con il fallback
        # e con hamming_within(); con e<=N sarebbero ammessi anche indel e
        # la posizione della CDR3 dedotta da m.start()+len(anchor) sarebbe sbagliata
        pattern = regex.compile(f"(?:{anchor}){{s<={max_mismatches}}}")

        def find_fuzzy(read_seq):
            m = pattern.search(read_seq)
            return m.start() if m else None
        return find_fuzzy

    L = len(anchor)

    def find_fallback(read_seq):
        limit = len(read_seq) - L
        for i in range(limit + 1):
            window = read_seq[i:i + L]
            mismatches = 0
            for a, b in zip(window, anchor):
                if a != b:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                return i
        return None
    return find_fallback


def hamming_within(seq_a, seq_b, max_mismatches):
    if len(seq_a) != len(seq_b):
        return False
    mismatches = 0
    for a, b in zip(seq_a, seq_b):
        if a != b:
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True


# ==========================================================================
# Traduzione
# ==========================================================================

def translate_cdr3(nt_seq):
    usable_len = (len(nt_seq) // 3) * 3
    if usable_len == 0:
        return None
    aa = str(Seq(nt_seq[:usable_len]).translate())
    if "*" in aa:
        return None
    return aa


# ==========================================================================
# Output CSV
# ==========================================================================

def write_counts_csv(path, counter, header):
    total = sum(counter.values())
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for seq, count in counter.most_common():
            freq = count / total if total else 0
            writer.writerow([seq, count, f"{freq:.6f}"])


def write_fasta_reads(path, records):
    """Scrive un FASTA con un record per ogni read valida (ridondante,
    cioe' senza deduplicare): utile per analisi che devono pesare per
    frequenza (es. logo plot con ggseqlogo)."""
    with open(path, "w") as f:
        for i, (read_id, seq) in enumerate(records, start=1):
            header = read_id if read_id else f"read_{i}"
            f.write(f">{header}\n{seq}\n")


# ==========================================================================
# Estrazione/conteggio CDR3 su un singolo fastq (gia' pronto/merged)
# ==========================================================================

def run_cdr3_extraction(fastq_path, ref_seq, args, out_prefix):
    """
    Esegue trimming + estrazione CDR3 basata su anchor + traduzione +
    conteggio su un singolo file fastq. Scrive i due CSV di output e
    ritorna un dizionario di statistiche riassuntive.
    """
    cdr3_start0 = args.cdr3_start - 1   # 0-based, inclusivo
    cdr3_end0 = args.cdr3_end           # 0-based, esclusivo
    cdr3_len = cdr3_end0 - cdr3_start0

    if not (0 <= cdr3_start0 < cdr3_end0 <= len(ref_seq)):
        sys.exit("Intervallo --cdr3-start/--cdr3-end non valido rispetto alla lunghezza del riferimento.")
    if cdr3_start0 < args.upstream_anchor_len:
        sys.exit("Il riferimento non ha abbastanza sequenza a monte per l'anchor richiesto "
                  "(riduci --upstream-anchor-len).")
    if cdr3_end0 + args.downstream_anchor_len > len(ref_seq):
        sys.exit("Il riferimento non ha abbastanza sequenza a valle per l'anchor richiesto "
                  "(riduci --downstream-anchor-len).")

    upstream_anchor = ref_seq[cdr3_start0 - args.upstream_anchor_len: cdr3_start0]
    downstream_anchor = ref_seq[cdr3_end0: cdr3_end0 + args.downstream_anchor_len]

    find_upstream = build_anchor_matcher(upstream_anchor, args.max_anchor_mismatches)
    trim_poly_g = build_poly_g_trimmer(0 if args.no_poly_g_trim else args.poly_g_min_run)

    n_reads_expected = count_fastq_reads(fastq_path)
    print(f"    [INFO] {n_reads_expected} reads da processare", file=sys.stderr)

    nt_counter = Counter()
    aa_counter = Counter()
    aa_reads = []  # (read_id, aa_seq) per ogni read con CDR3 tradotta senza stop codon

    stats = {
        "n_total": 0, "n_pass_trim": 0, "n_upstream_found": 0,
        "n_too_short": 0, "n_downstream_mismatch": 0,
        "n_cdr3_ok": 0, "n_translated_ok": 0,
    }

    with open_maybe_gzip(fastq_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            stats["n_total"] += 1

            if stats["n_total"] % 10000 == 0:
                pct = 100 * stats["n_total"] / n_reads_expected if n_reads_expected else 0
                print(f"\r    [INFO] {stats['n_total']}/{n_reads_expected} reads ({pct:.1f}%)",
                      end="", file=sys.stderr, flush=True)

            record = trim_poly_g(record)
            trimmed = sliding_window_trim(
                record, window=args.window, min_qual=args.min_qual, min_len=args.min_len)
            if trimmed is None:
                continue
            stats["n_pass_trim"] += 1

            read_seq = str(trimmed.seq).upper()

            up_pos = find_upstream(read_seq)
            if up_pos is None:
                continue
            stats["n_upstream_found"] += 1

            cdr3_start_in_read = up_pos + len(upstream_anchor)
            cdr3_end_in_read = cdr3_start_in_read + cdr3_len
            downstream_end_in_read = cdr3_end_in_read + len(downstream_anchor)

            if downstream_end_in_read > len(read_seq):
                stats["n_too_short"] += 1
                continue

            cdr3_nt = read_seq[cdr3_start_in_read:cdr3_end_in_read]
            downstream_seq = read_seq[cdr3_end_in_read:downstream_end_in_read]

            if not hamming_within(downstream_seq, downstream_anchor, args.max_anchor_mismatches):
                stats["n_downstream_mismatch"] += 1
                continue

            stats["n_cdr3_ok"] += 1
            nt_counter[cdr3_nt] += 1

            aa = translate_cdr3(cdr3_nt)
            if aa is None:
                continue
            stats["n_translated_ok"] += 1
            aa_counter[aa] += 1
            aa_reads.append((trimmed.id, aa))

    print("", file=sys.stderr)  # a capo dopo l'ultima riga di avanzamento

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_counts_csv(f"{out_prefix}_cdr3_nt_counts.csv", nt_counter,
                      header=("cdr3_nt", "count", "frequency"))
    write_counts_csv(f"{out_prefix}_cdr3_aa_counts.csv", aa_counter,
                      header=("cdr3_aa", "count", "frequency"))
    write_fasta_reads(f"{out_prefix}_cdr3_aa_reads.fasta", aa_reads)

    stats["n_unique_nt"] = len(nt_counter)
    stats["n_unique_aa"] = len(aa_counter)
    return stats


def print_stats(sample_name, stats):
    print(f"=== Riepilogo {sample_name} ===")
    print(f"  Reads totali:                    {stats['n_total']}")
    print(f"  Reads dopo quality trimming:      {stats['n_pass_trim']}")
    print(f"  Reads con anchor upstream:        {stats['n_upstream_found']}")
    print(f"  Scartate (troppo corte):          {stats['n_too_short']}")
    print(f"  Scartate (anchor downstream):     {stats['n_downstream_mismatch']}")
    print(f"  CDR3 estratte con successo:       {stats['n_cdr3_ok']}")
    print(f"  CDR3 tradotte senza stop codon:   {stats['n_translated_ok']}")
    print(f"  CDR3 nt uniche:                   {stats['n_unique_nt']}")
    print(f"  CDR3 aa uniche:                   {stats['n_unique_aa']}")


# ==========================================================================
# Orchestrazione batch: individuazione campioni, concatenazione lane, fastp
# ==========================================================================

SAMPLE_RE = re.compile(r"^(?P<sample>.+)_L\d+_R1_\d+\.fastq\.gz$")


def find_samples(raw_dir):
    raw_dir = Path(raw_dir)
    samples = set()
    for f in raw_dir.glob("*_R1_*.fastq.gz"):
        m = SAMPLE_RE.match(f.name)
        if m:
            samples.add(m.group("sample"))
    return sorted(samples)


def concatenate_lanes(lane_files, dest_path):
    """Concatena file .fastq.gz (multistream gzip valido) e verifica il conteggio reads."""
    with open(dest_path, "wb") as dest:
        for f in lane_files:
            with open(f, "rb") as src:
                shutil.copyfileobj(src, dest)

    n_sum = sum(count_fastq_reads(f) for f in lane_files)
    n_cat = count_fastq_reads(dest_path)
    if n_sum != n_cat:
        sys.exit(f"[ERRORE] QC concatenazione fallito per {dest_path}: "
                  f"somma lane={n_sum}, file concatenato={n_cat}")
    return n_cat


def ensure_fastp_available(auto_install=True):
    """
    Verifica che 'fastp' sia disponibile nel PATH. Se manca e auto_install e'
    True, prova a installarlo con 'conda install -c bioconda fastp'.
    """
    if shutil.which("fastp") is not None:
        return

    print("[INFO] 'fastp' non trovato nel PATH.", file=sys.stderr)

    if not auto_install:
        sys.exit("[ERRORE] 'fastp' non e' installato. Installalo con: "
                  "conda install -c bioconda fastp  (oppure usa --no-auto-install-fastp "
                  "solo se vuoi gestirlo manualmente e installarlo tu prima di rilanciare).")

    if shutil.which("conda") is None and shutil.which("mamba") is None:
        sys.exit("[ERRORE] 'fastp' non e' installato e ne' 'conda' ne' 'mamba' sono disponibili "
                  "per installarlo automaticamente. Installa fastp manualmente, ad es.: "
                  "conda install -c bioconda fastp")

    installer = "mamba" if shutil.which("mamba") is not None else "conda"
    print(f"[INFO] Provo a installare fastp con '{installer} install -c bioconda fastp'...",
          file=sys.stderr)
    try:
        subprocess.run([installer, "install", "-y", "-c", "bioconda", "fastp"], check=True)
    except subprocess.CalledProcessError:
        sys.exit(f"[ERRORE] Installazione automatica di fastp con {installer} fallita. "
                  f"Installalo manualmente con: conda install -c bioconda fastp")

    if shutil.which("fastp") is None:
        sys.exit("[ERRORE] fastp non risulta disponibile nel PATH anche dopo il tentativo "
                  "di installazione automatica (potrebbe servire riattivare l'ambiente conda).")

    print("[INFO] fastp installato correttamente.", file=sys.stderr)


def run_fastp_merge(cat_r1, cat_r2, merged_out, json_path, html_path, log_path, args):
    cmd = [
        "fastp",
        "-i", str(cat_r1), "-I", str(cat_r2),
        "--merge", "--merged_out", str(merged_out),
        "--overlap_len_require", str(args.overlap_len_require),
        "--overlap_diff_percent_limit", str(args.overlap_diff_percent_limit),
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--detect_adapter_for_pe",
        "--json", str(json_path),
        "--html", str(html_path),
    ]
    with open(log_path, "w") as log:
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        except FileNotFoundError:
            sys.exit("[ERRORE] 'fastp' non trovato nel PATH. Installalo (es. 'conda install -c "
                      "bioconda fastp') o usa la modalita' 'single' su un fastq gia' mergiato.")
        except subprocess.CalledProcessError:
            sys.exit(f"[ERRORE] fastp e' fallito per {cat_r1.name}/{cat_r2.name}. "
                      f"Controlla il log: {log_path}")


def run_batch(args):
    ensure_fastp_available(auto_install=not args.no_auto_install_fastp)

    raw_dir = Path(args.raw_dir)
    work_dir = Path(args.work_dir)
    results_dir = Path(args.results_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    ref_seq = load_reference(args.reference)

    samples = find_samples(raw_dir)
    if not samples:
        sys.exit(f"[ERRORE] Nessun campione trovato in {raw_dir} "
                  f"(atteso naming '<codice>_S<n>_L00X_R1_001.fastq.gz').")

    print(f"Trovati {len(samples)} campioni:")
    for s in samples:
        print(f"  - {s}")
    print()

    all_stats = []

    for sample in samples:
        print(f"=== Campione: {sample} ===")

        r1_files = sorted(raw_dir.glob(f"{sample}_L*_R1_*.fastq.gz"))
        r2_files = sorted(raw_dir.glob(f"{sample}_L*_R2_*.fastq.gz"))

        if len(r1_files) != args.n_lanes or len(r2_files) != args.n_lanes:
            print(f"[ATTENZIONE] {sample}: attese {args.n_lanes} lane R1/R2, "
                  f"trovate {len(r1_files)} R1 e {len(r2_files)} R2. Campione SALTATO.",
                  file=sys.stderr)
            continue

        cat_r1 = work_dir / f"{sample}_R1.fastq.gz"
        cat_r2 = work_dir / f"{sample}_R2.fastq.gz"
        merged = work_dir / f"{sample}_merged.fastq.gz"

        print(f"  -> Concateno {len(r1_files)} lane R1 e {len(r2_files)} lane R2")
        n_reads = concatenate_lanes(r1_files, cat_r1)
        concatenate_lanes(r2_files, cat_r2)
        print(f"  -> QC concatenazione OK: {n_reads} reads")

        print("  -> Merge R1+R2 con fastp")
        run_fastp_merge(
            cat_r1, cat_r2, merged,
            work_dir / f"{sample}_fastp.json",
            work_dir / f"{sample}_fastp.html",
            work_dir / f"{sample}_fastp.log",
            args,
        )

        print("  -> Estrazione/conteggio CDR3")
        out_prefix = results_dir / sample
        stats = run_cdr3_extraction(merged, ref_seq, args, out_prefix)
        print_stats(sample, stats)
        stats["sample"] = sample
        all_stats.append(stats)
        print()

    if all_stats:
        summary_path = results_dir / "summary_all_samples.csv"
        fieldnames = ["sample", "n_total", "n_pass_trim", "n_upstream_found",
                      "n_too_short", "n_downstream_mismatch", "n_cdr3_ok",
                      "n_translated_ok", "n_unique_nt", "n_unique_aa"]
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for s in all_stats:
                writer.writerow({k: s[k] for k in fieldnames})
        print(f"Riepilogo aggregato: {summary_path}")

    print(f"Tutti i campioni processati. Risultati in {results_dir}/")


def run_single(args):
    ref_seq = load_reference(args.reference)
    stats = run_cdr3_extraction(args.fastq, ref_seq, args, args.out_prefix)
    print_stats(Path(args.out_prefix).name, stats)
    print("File generati:")
    print(f"  {args.out_prefix}_cdr3_nt_counts.csv")
    print(f"  {args.out_prefix}_cdr3_aa_counts.csv")
    print(f"  {args.out_prefix}_cdr3_aa_reads.fasta")


# ==========================================================================
# CLI
# ==========================================================================

def add_cdr3_analysis_args(parser):
    parser.add_argument("--reference", required=True,
                         help="FASTA con la sequenza 'default' pIII-scFv (nucleotidica, in frame)")
    parser.add_argument("--cdr3-start", type=int, default=361,
                         help="Posizione iniziale CDR3 nel riferimento, 1-based, in nucleotidi (default 361)")
    parser.add_argument("--cdr3-end", type=int, default=402,
                         help="Posizione finale CDR3 nel riferimento, 1-based, inclusiva, in nucleotidi (default 402)")
    parser.add_argument("--upstream-anchor-len", type=int, default=20,
                         help="Lunghezza (nt) dell'anchor costante a monte della CDR3 (default 20)")
    parser.add_argument("--downstream-anchor-len", type=int, default=20,
                         help="Lunghezza (nt) dell'anchor costante a valle della CDR3 (default 20)")
    parser.add_argument("--max-anchor-mismatches", type=int, default=2,
                         help="Mismatch (solo sostituzioni) massime tollerate per ciascun anchor, "
                              "upstream e downstream (default 2; 0 = match esatto)")
    parser.add_argument("--min-qual", type=int, default=20,
                         help="Qualita' media minima per la finestra di trimming (default 20)")
    parser.add_argument("--window", type=int, default=4,
                         help="Ampiezza finestra scorrevole per il trimming (default 4)")
    parser.add_argument("--min-len", type=int, default=50,
                         help="Lunghezza minima della read dopo il trimming (default 50)")
    parser.add_argument("--poly-g-min-run", type=int, default=10,
                         help="Lunghezza minima di una coda di poli-G da rimuovere prima del "
                              "quality trimming (default 10, tipico di NextSeq/NovaSeq)")
    parser.add_argument("--no-poly-g-trim", action="store_true",
                         help="Disattiva la rimozione delle code di poli-G")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    p_batch = subparsers.add_parser(
        "batch", help="Processa tutti i campioni multi-lane trovati in una cartella")
    p_batch.add_argument("--raw-dir", required=True,
                          help="Cartella con i fastq.gz grezzi di tutti i campioni")
    p_batch.add_argument("--work-dir", default="work",
                          help="Cartella per i file intermedi (concatenati, merged) (default: work)")
    p_batch.add_argument("--results-dir", default="results",
                          help="Cartella per i CSV finali per campione (default: results)")
    p_batch.add_argument("--n-lanes", type=int, default=4,
                          help="Numero di lane attese per campione (default 4)")
    p_batch.add_argument("--overlap-len-require", type=int, default=30,
                          help="fastp: overlap minimo richiesto per il merge R1+R2 (default 30)")
    p_batch.add_argument("--overlap-diff-percent-limit", type=int, default=10,
                          help="fastp: percentuale massima di mismatch tollerata nell'overlap (default 10)")
    p_batch.add_argument("--no-auto-install-fastp", action="store_true",
                          help="Non tentare di installare fastp automaticamente con conda/mamba "
                               "se manca nel PATH (fallisce e basta con un messaggio d'errore)")
    add_cdr3_analysis_args(p_batch)
    p_batch.set_defaults(func=run_batch)

    p_single = subparsers.add_parser(
        "single", help="Processa un singolo fastq gia' pronto (es. gia' mergiato)")
    p_single.add_argument("--fastq", required=True, help="FASTQ da analizzare (anche .gz)")
    p_single.add_argument("--out-prefix", required=True, help="Prefisso dei file di output")
    add_cdr3_analysis_args(p_single)
    p_single.set_defaults(func=run_single)

    args = parser.parse_args()

    if args.max_anchor_mismatches < 0:
        parser.error("--max-anchor-mismatches deve essere >= 0")

    if not HAVE_REGEX and args.max_anchor_mismatches > 0:
        print("[WARN] modulo 'regex' non installato: uso il fallback Python puro per la "
              "ricerca approssimata degli anchor (piu' lento su dataset grandi). "
              "Per velocizzare: pip install regex", file=sys.stderr)

    args.func(args)


if __name__ == "__main__":
    main()
