# Scatter plot log-log tra due analisi, con retta 1:1,
# regressione di potenza (y = a * x^b) e R^2
# ============================================================

library(readxl)
library(ggplot2)
library(scales)
library(stringr)

# ---------- 1. PARAMETRI DA ADATTARE ----------
file_input <- "~/input pII-pIII comparison.xlsx"      # percorso del file
col_x      <- "pII-pIII enrichment"       # nome colonna x
col_y      <- "MM473-MET5a enrichment"       # nome colonna y
file_png   <- "grafico.png"    # nome del grafico salvato

# ---------- 2. LETTURA ----------
dati <- read_excel(file_input)
# Per un CSV con virgola decimale usa invece:
# dati <- read.csv2("file.csv")

df <- data.frame(
  x = dati[[col_x]],
  y = dati[[col_y]]
)

# ---------- 3. PULIZIA ----------
# Controllo delle celle non numeriche (es. "unique") prima della conversione
non_num_x <- df$x[is.na(suppressWarnings(as.numeric(gsub(",", ".", df$x)))) & !is.na(df$x)]
non_num_y <- df$y[is.na(suppressWarnings(as.numeric(gsub(",", ".", df$y)))) & !is.na(df$y)]
cat("Valori non numerici in x:\n"); print(table(non_num_x))
cat("Valori non numerici in y:\n"); print(table(non_num_y))

# Conversione in numerico ("unique" e simili diventano NA)
df$x <- suppressWarnings(as.numeric(gsub(",", ".", df$x)))
df$y <- suppressWarnings(as.numeric(gsub(",", ".", df$y)))

n_tot <- nrow(df)

# Teniamo solo le sequenze presenti in entrambe le analisi (valori > 0)
df <- df[!is.na(df$x) & !is.na(df$y) & df$x > 0 & df$y > 0, ]
cat(sprintf("Sequenze totali: %d | usate nel grafico e nella regressione: %d\n",
            n_tot, nrow(df)))

# ---------- 4. REGRESSIONE IN SCALA LOG-LOG ----------
fit <- lm(log10(y) ~ log10(x), data = df)
a   <- unname(coef(fit)[1])
b   <- unname(coef(fit)[2])
r2  <- summary(fit)$r.squared

cat(sprintf("y = %.4g * x^%.4f | R2 (scala log) = %.4f\n", 10^a, b, r2))

etichetta <- sprintf("y = %.3g · x^%.3f\nR² = %.3f\nn = %d",
                     10^a, b, r2, nrow(df))

# ---------- 5. GRAFICO ----------
numeri_chiari <- function(v) format(v, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)

titolo <- "Scattering plot of the common sequences between the third pannings of MET5\u03b1 and MM473 against the enrichment rate from the previous panning displays a weak positive correlation"

p <- ggplot(df, aes(x, y)) + 
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "blue", se = FALSE) +
  scale_x_log10(labels = numeri_chiari, breaks = log_breaks(n = 8)) +
  scale_y_log10(labels = numeri_chiari, breaks = log_breaks(n = 8)) +
  annotate("text", x = min(df$x), y = max(df$y), label = etichetta,
           hjust = 0, vjust = 1, size = 4) +
  labs(
    title = str_wrap(titolo, width = 70),
    x = "MM473 pII-pIII enrichment",
    y = "MM473-MET5\u03b1 enrichment"
  ) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.title.position = "plot"
  )
print(p)
