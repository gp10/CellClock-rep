Plot_sbEdU_tEdU_dotplot <- function() {
  ggplot(mydata, aes(x=factor(Genotype, levels=c("WT","Mut")), y=sbEdU_tEdU, color = Animal)) +
    geom_jitter(aes(shape = Animal), color = "darkgray", position = position_jitter(0.2),size = 1.2) +
    scale_shape_manual(values=c(0:6)) + 
    geom_dotplot(data = mydata_avg,aes(shape = Animal), color = "black",
                 binaxis='y', stackdir='center', dotsize = 0.8, position = position_dodge(0.3)) +
    stat_summary(data = mydata_avg, fun.data = "mean_sdl", fun.args = list(mult=1),
                 geom = "pointrange", color = "red") +
    ylab("sb EdU / total EdU") + xlab("Genotype") + ylim(0.3,0.5) 
  #a + geom_jitter(data = mydata_avg, aes(fill = Animal, shape = Animal), color = "black", position = position_jitter(0.2),size = 2) +
  #  scale_fill_manual(values=c(rep("black",7)))
}