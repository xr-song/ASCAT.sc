getStartsEnds <- function (window, chr, lengthChr, dna = NULL, pathWindows = NA, 
    pathBadBins = NA, centromeres = NULL, excludebed = NULL) 
{
    if (!is.na(pathWindows)) {
        badbins <- if (!is.na(pathBadBins)) 
            read.table(pathBadBins)[, 1]
        else NULL
        t <- read.table(pathWindows, header = TRUE)
        subt <- t[t$CHR == chr, ]
        starts <- c(1, as.numeric(as.character(subt$END[-length(subt$END)])) + 1)
        ends <- as.numeric(as.character(subt$END))
        # modification: exclude intervals after starts/ends are computed to avoid wrong interval in the first position of chromosome
        if (!is.null(badbins) && length(badbins) > 0) {
            # badbins refers to row indices in original t, but we subset to subt
            # so need to map badbins to indices within subt, if necessary
            # assuming badbins are relative to the full t, restrict to this chr
            chr_row_indices <- which(t$CHR == chr)
            badbins_in_chr <- which(chr_row_indices %in% badbins)
            if(length(badbins_in_chr) > 0) {
                starts <- starts[-badbins_in_chr]
                ends <- ends[-badbins_in_chr]
            }
        }
    }
    else if (!is.null(excludebed)) {
        keepChr <- gsub("chr", "", excludebed[, 1]) == gsub("chr", "", chr)
        if (sum(keepChr) > 0) {
            excludebed <- excludebed[keepChr, , drop = FALSE]
            excludebed <- rbind(excludebed, excludebed[nrow(excludebed), ])
            excludebed[nrow(excludebed), 2:3] <- lengthChr - 1:0
            grE <- GRanges(chr, IRanges(excludebed[, 2], excludebed[, 3]))
            covs <- coverage(grE)
            effectiveLength <- sum(covs == 0)
            divideChr <- seq(1, effectiveLength, window)
            starts <- divideChr[-length(divideChr)]
            effectivePos <- which(as.logical((covs == 0)[[1]]))
            starts <- effectivePos[starts]
            ends <- c(starts[-1] - 1, lengthChr)
        }
        else {
            return(getStartsEnds(window = window, chr = chr, 
                lengthChr = lengthChr, centromeres = centromeres))
        }
    }
    else {
        divideChr <- seq(0, lengthChr, window)
        starts <- divideChr[-length(divideChr)] + 1
        ends <- divideChr[-1]
    }

    excludeBadBins(
        removeCentromeres(list(starts = starts, ends = ends), chr = chr, centromeres = centromeres), 
        chr = gsub("chr", "", chr), dna = dna
    )
}



writeLSeToBed <- function(lSe, file) {
  bed_list <- list()

  for (chrom in names(lSe)) {
    if (chrom %in% c("starts", "ends")) next 
    starts <- lSe[[chrom]]$starts
    ends <- lSe[[chrom]]$ends
    if (is.null(starts) || is.null(ends)) {
      starts <- lSe[[paste0("starts")]]
      ends <- lSe[[paste0("ends")]]
      if (is.null(starts) || is.null(ends)) next
    }
    bed_df <- data.frame(
      chr = chrom,
      start = as.numeric(starts),
      end = as.numeric(ends)
    )
    bed_list[[chrom]] <- bed_df
  }

  bed_df_all <- do.call(rbind, bed_list)
  write.table(
    bed_df_all, file = file, sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
    message(paste('Written to', file))
}
