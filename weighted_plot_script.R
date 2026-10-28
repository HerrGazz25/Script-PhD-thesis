# ============================================================
# Sequence logo pesato per read count - Top 30 CDR3 arricchite
# ============================================================
# Costruisce un logo che riflette l'abbondanza reale (readCount)
# di ciascuna CDR3, aggirando il limite di Seq2Logo (che tratta
# l'header FASTA come semplice etichetta e non pesa per il conteggio).
#
# Header FASTA atteso: >qNp7PzQZ1ASQtwclTE04e|readCount_31417

# ---- 1. Pacchetti richiesti -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("ggseqlogo", quietly = TRUE)) install.packages("ggseqlogo")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(Biostrings)
library(ggseqlogo)
library(stringr)
library(ggplot2)

# ---- 2. Import del FASTA -----------------------------------
# Sostituisci il path con quello del tuo file
fasta_path <- "pr1.fasta"
seqs <- readAAStringSet(fasta_path)

headers  <- names(seqs)
aa_seqs  <- as.character(seqs)

# ---- 3. Estrazione del read count dall'header --------------
# Header nel formato: >ID|readCount_31417
read_counts <- as.numeric(str_extract(headers, "(?<=readCount_)\\d+"))

if (any(is.na(read_counts))) {
  stop("Impossibile estrarre readCount da uno o più header: controlla il formato.")
}

# ---- 4. Controllo di qualità sulle sequenze -----------------
cdr3_length <- 14L
if (!all(nchar(aa_seqs) == cdr3_length)) {
  stop("Trovate sequenze con lunghezza diversa da ", cdr3_length,
       " aa: verifica l'allineamento prima di procedere.")
}

# ---- 5. Costruzione della matrice di frequenza pesata (PFM) --
aa_alphabet <- c("A","R","N","D","C","Q","E","G","H","I",
                 "L","K","M","F","P","S","T","W","Y","V")

pfm <- matrix(0, nrow = length(aa_alphabet), ncol = cdr3_length,
              dimnames = list(aa_alphabet, seq_len(cdr3_length)))

for (i in seq_along(aa_seqs)) {
  residues <- strsplit(aa_seqs[i], "")[[1]]
  w <- read_counts[i]
  for (pos in seq_len(cdr3_length)) {
    aa <- residues[pos]
    if (aa %in% aa_alphabet) {
      pfm[aa, pos] <- pfm[aa, pos] + w
    }
  }
}

# Facoltativo: comprimere il range di abbondanza (da 2 a 31417)
# per evitare che 1-2 cloni dominino totalmente il logo.
# Decommenta la riga seguente per usare il log10(readCount+1) come peso
# al posto del readCount grezzo:
# pfm_log <- matrix(0, nrow = length(aa_alphabet), ncol = cdr3_length,
#                    dimnames = list(aa_alphabet, seq_len(cdr3_length)))
# for (i in seq_along(aa_seqs)) {
#   residues <- strsplit(aa_seqs[i], "")[[1]]
#   w <- log10(read_counts[i] + 1)
#   for (pos in seq_len(cdr3_length)) {
#     aa <- residues[pos]
#     if (aa %in% aa_alphabet) pfm_log[aa, pos] <- pfm_log[aa, pos] + w
#   }
# }

# ---- 6. Generazione del logo --------------------------------
# method = "bits"        -> altezza totale = information content (0-4.32 bit)
# method = "probability" -> altezza totale sempre 1, mostra solo composizione
p <- ggseqlogo(pfm, method = "bits", seq_type = "aa") +
  xlab("Position CDR3 loop") +
  ylab("Bits") +
  theme_logo(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

print(p)

# ---- 7. Esportazione ad alta risoluzione ---------------------
ggsave("weighted_s1.png", plot = p, width = 10, height = 4, dpi = 600)
ggsave("CDR3_logo_weighted.pdf", plot = p, width = 10, height = 4)

