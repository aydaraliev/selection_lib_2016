library(GenWin)

extract_snps <- function(computed_windows, chr) {
  snps_list = list()
  for (i in 1:nrow(computed_windows)) {
    row = computed_windows[i,]
    end_coordinate = row$WindowStop
    start_coordinate = row$WindowStart
    extracted_snps = subset(chr, chr$Position < end_coordinate & chr$Position > start_coordinate)
    snps_list[[i]] = extracted_snps
  }
  snps_list = do.call(rbind, snps_list)
  return(snps_list)
}
  
chr_top_snps_list = list()

for (chr in 1:22) {
  chr = fst[fst$Chromosome == chr,]
  chr = chr[order(chr$Position),]
  chr = chr[!duplicated(chr$Position),]
  #chr$Overall.Fst = scale(chr$Overall.Fst)
  raw_fst_quant = quantile(chr$Overall.Fst, 0.999)
  chrom_num = chr$Chromosome[1]
  jpeg(filename = paste(chrom_num, ".jpg", sep=""), width = 1200, height = 800)
  chr_comp = splineAnalyze(chr$Overall.Fst, chr$Position, plotRaw = 1, plotWindows = 1, method = 4)
  wstat_vector = chr_comp$windowData$Wstat[!is.na(chr_comp$windowData$Wstat)]
  abline(quantile(wstat_vector, 0.999), 0, col = "green")
  title(main = paste(round(quantile(wstat_vector, 0.999), 3), "threshold (0.999 quantile)", sep = " "))
  window_q = quantile(wstat_vector, 0.999)
  chr_wstat_no_na = chr_comp$windowData[!is.na(chr_comp$windowData$Wstat),]
  chr_top_windows = chr_wstat_no_na[chr_wstat_no_na$Wstat > window_q,]
  chr_top_snps_list[[chrom_num]] = extract_snps(chr_top_windows, chr)
  layout(matrix(c(2,1), ncol = 1))
  par(new = T)
  plot(Overall.Fst ~ Position, chr, ann = F, type = "n", axes = F)
  first_graph_threshold = quantile(chr$Overall.Fst[!is.na(chr$Overall.Fst)], 0.999)
  abline(first_graph_threshold, 0, col = "green")
  title(main = paste(round(first_graph_threshold, 3), "threshold (0.999 quantile)", sep = " "))
  write.table(chr_comp$windowData, file = paste(chrom_num, ".txt", sep=""), sep = '\t')
  dev.off()
}

chr_top_snps_list = do.call(rbind, chr_top_snps_list)
write.csv(chr_top_snps_list, file = "top_snps.txt", sep = '\t')

