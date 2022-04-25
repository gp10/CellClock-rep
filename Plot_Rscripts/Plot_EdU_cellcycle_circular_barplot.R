Plot_EdU_cellcycle_circular_barplot <- function() {
  # library
  library(tidyverse)
  
  # Import dataset (done)
  # freq_bEdU_avg
  # Order data (optional):
  # freq_bEdU_avg = freq_bEdU_avg %>% arrange(Phase, Fraction)
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 2
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(freq_bEdU_avg$Phase), ncol(freq_bEdU_avg)) )
  colnames(to_add) <- colnames(freq_bEdU_avg)
  to_add$Phase <- rep(levels(freq_bEdU_avg$Phase), each=empty_bar)
  freq_bEdU_avg <- rbind(freq_bEdU_avg, to_add)
  freq_bEdU_avg <- freq_bEdU_avg %>% arrange(Phase)
  freq_bEdU_avg$id <- seq(1, nrow(freq_bEdU_avg))
  
  # Get the name and the y position of each label
  # label_data <- data
  # number_of_bar <- nrow(label_data)
  # angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  # label_data$hjust <- ifelse( angle < -90, 1, 0)
  # label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data <- freq_bEdU_avg %>% 
    group_by(Phase) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  #grid_data <- grid_data[-1,]
  grid_data$start[1] = grid_data$end[1]+empty_bar-1
  
  # Make the plot
  #pdf("EdU_cellcycle_circular_barplot.pdf")
  p <- ggplot(freq_bEdU_avg, aes(x=as.factor(id), y=Fraction, fill=Genotype)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=Fraction, fill=Genotype), stat="identity", alpha=0.5) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(freq_bEdU_avg$id),6), y = c(0, 20, 40, 60, 80, 100), label = c("0","20", "40", "60", "80", "100") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
    
    geom_bar(aes(x=as.factor(id), y=Fraction, fill=Genotype), stat="identity", alpha=0.5) +
    scale_fill_manual("legend", values = c("WT" = "black", "Mut" = "red")) +
    
    ylim(-100,120) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() + 
    #geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=c("G1", "S", "G2", "M")), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
  
  p
  #dev.off()
}
