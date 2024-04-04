library(VariantAnnotation)
library(tidyverse)
library(plyranges)
f <- "vcf/1000GP_SNV_INDEL_panhuman_unrel_0.5.vcf.gz"
vcf <- VcfFile(f)
sq <- seqinfo(vcf)
isCircular(sq) <- rep(FALSE, length(seqnames(sq)))
genome(sq) <- "CRCh38"
## Try refining the import to be viable
param <- ScanVcfParam(
    fixed = c("ALT"), info = NA
    # which = GRanges(sq, seqinfo = sq)[21:22,]
)
gr <- readVcf(vcf, genome = sq, param = param) |>
    rowRanges() |>
    select(REF, ALT)
## Can we quickly categorise the variant types
gr %>%
    mutate(
        alt_len = unlist(nchar(ALT)),
        type = case_when(
            width == 1 & alt_len == 1 ~ "SNV",
            alt_len > width ~ "Ins",
            alt_len < width ~ "Del",
        )
    ) %>%
    mcols() %>%
    .[["type"]] %>%
    table()

x <- as.character(unlist(gr$ALT)[1:100])
unlist(nchar(x))
