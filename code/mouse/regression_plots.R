library(ggplot2)


topGeneList <- read.csv("../FinalData/wgr_annotated_top_rank_genelist.csv")

ggplot(topGeneList, aes(x=topGeneList$LHGV, y=topGeneList$WGR.Score)) +
  geom_point(shape=1)+
  stat_smooth(method="rlm")

library(slidify)
author("Lhgv1.33")
