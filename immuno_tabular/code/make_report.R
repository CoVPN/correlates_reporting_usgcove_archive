library(rmarkdown)
library(here)

rmd <- "report.Rmd"
outfile <- sprintf("%s-%s>", rmd, Sys.Date())

# Render pdf
render(
  input = rmd,
  output_format = "all", # Generate all outputs in YAML header
  output_dir = here("output")
)


