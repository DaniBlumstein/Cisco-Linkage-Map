#install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
library(ggthemes)
jpeg("male_map.jpg",height=8,width=10,units="in",res=700)
#Customize plot
par(fig=c(0,1,0,1),mar=c(1,1,1,1))

LinkageMap<-read.delim("F:/LepMap3/diploids_for_lepmap/final_male_map/diploid_combo_final_lepmap_separate5_order_R.txt")

#head(LinkageMap)
map <- ggplot(LinkageMap, aes(x = LG_NEW, y = male_position)) +
  geom_linerange(aes(x = LG_NEW, ymin = 0, ymax = male_position), inherit.aes = FALSE) + 
  geom_point(size =2, alpha=.5, color = "black")  +
  labs(y = "cM", x = "Linkage Groups", colour="") +
  ggtitle("Male Linkage Map") +
  scale_x_continuous(breaks = seq(1,40,3)) +
  scale_y_continuous(breaks = seq(0,200,10)) +
  theme_base() +
  theme(axis.text = element_text(size=16), axis.title= element_text(size= 18))+
  theme(legend.position = "bottom", legend.text = element_text(size=13))+
  theme(plot.title = element_text(hjust = 0.5, size= 18, face="bold"))
map

dev.off()




