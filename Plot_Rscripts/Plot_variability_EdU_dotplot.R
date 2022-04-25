# Dot plot colored by group
Plot_variability_EdU_dotplot <- function() {
  par(mfrow=c(2,2))    # set the plotting area into a 1*2 array
  require(gridExtra)
  plot1 <- ggplot(mydata, aes(x=Genotype, y=(EdU_blank/bEdU), fill=Animal)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.4), color="red") +
    ylim(0,1)
  plot1 <- plot1 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                geom="pointrange", color="black", size=0.1, position=position_dodge(0.4))
  
  plot2 <- ggplot(mydata, aes(x=Genotype, y=(EdU_BrdU/bEdU), fill=Animal)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.4), color="red") +
    ylim(0,1)
  plot2 <- plot2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                geom="pointrange", color="black", size=0.1, position=position_dodge(0.4))
  
  plot3 <- ggplot(mydata, aes(x=Genotype, y=(EdU_pH3G2/bEdU), fill=Animal)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.4), color="red") +
    ylim(0,1)
  plot3 <- plot3 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                geom="pointrange", color="black", size=0.1, position=position_dodge(0.4))
  
  plot4 <- ggplot(mydata, aes(x=Genotype, y=(EdU_pH3M/bEdU), fill=Animal)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.4), color="red") +
    ylim(0,1)
  plot4 <- plot4 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                geom="pointrange", color="black", size=0.1, position=position_dodge(0.4))
  grid.arrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2)
}