Plot_cellcycle_errorbars <- function() {
  # Variability_errorbar plot:
  freq_bAllandbEdU_avgmice$Phase <- factor(freq_bAllandbEdU_avgmice$Phase,levels=c("pH3M","pH3G2","BrdU","Blank"))
  freq_bAllandbEdU_avgmice$up[freq_bAllandbEdU_avgmice$Phase == "BrdU"] <- with(freq_bAllandbEdU_avgmice,
                                                                                up[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  freq_bAllandbEdU_avgmice$dn[freq_bAllandbEdU_avgmice$Phase == "BrdU"] <- with(freq_bAllandbEdU_avgmice,
                                                                                dn[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  freq_bAllandbEdU_avgmice$up[freq_bAllandbEdU_avgmice$Phase == "pH3G2"] <- with(freq_bAllandbEdU_avgmice,
                                                                                 up[Phase == "pH3G2"] + Fraction[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  freq_bAllandbEdU_avgmice$dn[freq_bAllandbEdU_avgmice$Phase == "pH3G2"] <- with(freq_bAllandbEdU_avgmice,
                                                                                 dn[Phase == "pH3G2"] + Fraction[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  freq_bAllandbEdU_avgmice$up[freq_bAllandbEdU_avgmice$Phase == "pH3M"] <- with(freq_bAllandbEdU_avgmice,
                                                                                up[Phase == "pH3M"] + Fraction[Phase == "pH3G2"] + Fraction[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  freq_bAllandbEdU_avgmice$dn[freq_bAllandbEdU_avgmice$Phase == "pH3M"] <- with(freq_bAllandbEdU_avgmice,
                                                                                dn[Phase == "pH3M"] + Fraction[Phase == "pH3G2"] + Fraction[Phase == "BrdU"] + Fraction[Phase == "Blank"])
  
  dodge <- position_dodge(width=-0.25)
  ggplot(freq_bAllandbEdU_avgmice, aes(x=factor(Lineage), y=Fraction, fill=Phase)) +
    facet_grid(. ~ factor(Genotype, levels = c("WT","Mut"))) +
    theme(legend.position="right") +
    geom_bar(position="stack", stat="identity") +
    geom_errorbar(aes(ymin=dn, ymax=up),
                  width=0.5,                    # Width of the error bars
                  #position = position_dodge2(width = 0.5, padding = 0.5))
                  position=dodge)
}