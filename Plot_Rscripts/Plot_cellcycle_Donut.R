Plot_cellcycle_Donut <- function() {
  # load library
  library(ggplot2)
  library(patchwork)
  
  # Create test data.
  freq_bEdU_avgmice_mut <- freq_bEdU_avgmice[freq_bEdU_avgmice$Genotype=="Mut",]
  freq_bEdU_avgmice_wt <- freq_bEdU_avgmice[freq_bEdU_avgmice$Genotype=="WT",]
  
  # Compute percentages
  #data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  freq_bEdU_avgmice_mut$ymax <- cumsum(freq_bEdU_avgmice_mut$Fraction)
  freq_bEdU_avgmice_wt$ymax <- cumsum(freq_bEdU_avgmice_wt$Fraction)
  
  # Compute the bottom of each rectangle
  freq_bEdU_avgmice_mut$ymin <- c(0, head(freq_bEdU_avgmice_mut$ymax, n=-1))
  freq_bEdU_avgmice_wt$ymin <- c(0, head(freq_bEdU_avgmice_wt$ymax, n=-1))
  
  # Compute label position
  freq_bEdU_avgmice_mut$labelPosition <- (freq_bEdU_avgmice_mut$ymax + freq_bEdU_avgmice_mut$ymin) / 2
  freq_bEdU_avgmice_wt$labelPosition <- (freq_bEdU_avgmice_wt$ymax + freq_bEdU_avgmice_wt$ymin) / 2
  
  # Compute a good label
  #freq_bEdU_avgmice_mut$label <- paste0(freq_bEdU_avgmice_mut$Phase, "\n value: ", freq_bEdU_avgmice_mut$count)
  
  # Make the plot
  p <- ggplot(freq_bEdU_avgmice_mut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phase)) +
    geom_rect(color="white") +
    ggtitle("Mut") +
    #geom_text( x=2, aes(y=labelPosition, label=label, color=Phase), size=6) + # x here controls label position (inner / outer)
    #scale_fill_brewer(palette=3) +
    #scale_color_brewer(palette=3) +
    coord_polar(theta="y") +
    xlim(c(0, 4)) + # If xlim left boundary is big, no empty circle. If xlim is low, the ring becomes thinner.
    theme_void()
  #theme(legend.position = "none")
  q <- ggplot(freq_bEdU_avgmice_wt, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phase)) +
    geom_rect(color="white") +
    ggtitle("WT") +
    #geom_text( x=2, aes(y=labelPosition, label=label, color=Phase), size=6) + # x here controls label position (inner / outer)
    #scale_fill_brewer(palette=3) +
    #scale_color_brewer(palette=3) +
    coord_polar(theta="y") +
    xlim(c(0, 4)) + # If xlim left boundary is big, no empty circle. If xlim is low, the ring becomes thinner.
    theme_void() +
    theme(legend.position = "none")
  
  q + p
}