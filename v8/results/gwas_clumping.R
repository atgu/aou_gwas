library(data.table)
library(parallel)

prefix <- "gs://aou_wlu/v8_analysis/clumped_results"
pops <- c("AFR","AMR","EAS","EUR", "MID", "SAS", "META")  # <-- edit to your set

parse_meta <- function(uri) {
  pop <- sub("^.*?/clumped_results/([^/]+)/.*$", "\\1", uri)
  base <- basename(uri)
  chrom <- sub("^.*_chr([0-9]+|X|Y|M)\\.clumped$", "\\1", base)
  pheno <- sub(paste0("^", pop, "_(.*)_chr([0-9]+|X|Y|M)\\.clumped$"), "\\1", base)
  list(POP = pop, PHENO = pheno, CHROM = paste0("chr", chrom))
}

read_one <- function(uri, i, n) {
  message(sprintf("  [%d/%d] Reading %s", i, n, basename(uri)))
  meta <- parse_meta(uri)
  dt <- fread(sprintf("gcloud storage cat %s", shQuote(uri)))
  dt[, `:=`(ancestry = meta$POP, Phenotype = meta$PHENO, chromosome = meta$CHROM)]
  setcolorder(dt, c("ancestry","Phenotype","chromosome",
                    setdiff(names(dt), c("ancestry","Phenotype","chromosome"))))
  dt
}

all <- list()

for (pop in pops) {
  message(sprintf("Listing files for POP = %s ...", pop))
  uris <- system(
    sprintf("gcloud storage ls %s/%s/results/phenotype_*/*.clumped", prefix, pop),
    intern = TRUE
  )

  if (!length(uris)) {
    message(sprintf("  No files found for %s", pop))
    next
  }

  message(sprintf("  Found %d files for %s", length(uris), pop))

  all[[pop]] <- rbindlist(
    lapply(seq_along(uris), function(i) read_one(uris[i], i, length(uris))),
    fill = TRUE
  )
}

dt <- rbindlist(all, fill = TRUE)
fwrite(dt, "~/Desktop/all_clumped_annotated.tsv", sep = "\t")
