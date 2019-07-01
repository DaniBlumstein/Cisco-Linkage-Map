#install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
library(ggthemes)
jpeg("female_map.jpg",height=8,width=10,units="in",res=700)
#Customize plot
par(fig=c(0,1,0,1),mar=c(1,1,1,1))

centro <- read.delim("F:/LepMap3/consensus_files_for_lepmap/final_lepmap/centro_meta_only.txt")

LinkageMap<-read.delim("F:/LepMap3/consensus_files_for_lepmap/final_lepmap/consensus_hap_dip_lepmap3_lod15_final_renamed.txt")

map <- ggplot(LinkageMap, aes(x = LG_NEW, y = female_position, color=DuplicationStatus), centro) +
  geom_point(data = centro, aes(x =group, y=pos, size = 100, color="Centromere")) +
  geom_linerange(aes(x = LG_NEW, ymin = 0, ymax = female_position), inherit.aes = FALSE) + 
  geom_point(size =2, alpha=.5)  +
  scale_color_manual(values = c("red", "blue", "grey 60" )) +
  ggtitle("Female Linkage Map:" + "/n" + "Total Loci: 20,458, Duplicated Loci: 3,390, Singleton Loci: 17,068") +
  labs(y = "cM", x = "Linkage Groups", colour="") +
  scale_x_continuous(breaks = seq(1,38,2)) +
  scale_y_continuous(breaks = seq(0,200,10)) +
  theme_base() +
  theme(axis.text = element_text(size=16), axis.title= element_text(size= 18))+
  theme(legend.position = "top", legend.text = element_text(size=13))+
  theme(plot.title = element_text(hjust = 0.5, size= 18, face="bold")) +
  guides(size=FALSE)
map 

dev.off()
