# 08: render paper/paper.qmd (IEEE format) to PDF and HTML with Quarto.
# The conda environment ships a non-working quarto shim, so the standalone
# Quarto 1.5.57 installation is preferred when present. A UTF-8 locale is set
# so that symbols such as the minus sign and micro sign survive in tables.
Sys.setenv(LANG = "C.UTF-8", LC_ALL = "C.UTF-8")
quarto <- if (file.exists("/opt/quarto/opt/quarto/bin/quarto")) "/opt/quarto/opt/quarto/bin/quarto" else Sys.which("quarto")
for (fmt in c("pdf", "html")) {
  st <- system2(quarto, c("render", "paper/paper.qmd", "--to", fmt))
  if (st != 0) warning("Quarto render to ", fmt, " failed")
}
