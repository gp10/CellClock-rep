Plot_cellcycle_Piechart <- function() {
  # Prepare raw count table for pie chart:
  freq_bEdU_avgmice_count <- freq_bEdU_avgmice[rep(seq_len(nrow(freq_bEdU_avgmice)),
                                                   times=freq_bEdU_avgmice$n_FoV),
                                               -ncol(freq_bEdU_avgmice)]
    
  ### PLOT Pie Charts:
  p <- ggplot(data=freq_bEdU_avgmice_count, aes(x=factor(1), stat="bin", fill=Phase)) +
    theme_void() +
    geom_bar(position="fill", color="white") # Stacked bar chart
  p <- p + ggtitle("Phase by genotype") + xlab("") + ylab("Phase") # Adds titles
  p <- p + facet_grid(facets= . ~ Genotype) # Side by side bar chart
  p <- p + coord_polar(theta="y", start=0, direction=-1) # side by side pie chart
  p
}