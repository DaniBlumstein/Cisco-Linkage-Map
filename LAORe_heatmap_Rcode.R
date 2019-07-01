library(readxl)
library(ggplot2)
library(magrittr)
library(tidyr)

install.packages("viridis")
library(viridis)

jpeg("heatmap.jpg",height=8,width=10,units="in",res=150)
homeologs<-read_excel("Copy of LAORe_table_v2.xlsx", sheet="flipped table")
#convert to long table
homeologs_long<-homeologs %>% gather(Chromosome, Rank, -Species, -avg)
#rename columns
colnames(homeologs_long)<-c("Proto_Karyotype","avg","Species","Rank")
#convert to long file  
homeologs_long<-homeologs %>% gather(Species, Rank, -Proto_Karyotype, -avg)
#set rank as numeric for color coding
homeologs_long$Rank<-as.numeric(homeologs_long$Rank)
#convert species to factor for ordering plot axis, use same order as excel file
homeologs_long$Species<-factor(homeologs_long$Species, levels=colnames(homeologs)[2:13])
#convert protokaryotypes to factor for ordering plot axis, use same order as excel file
homeologs_long$Proto_Karyotype<-factor(homeologs_long$Proto_Karyotype, levels=rev(homeologs$Species))

#plot using geom_raster
homeologs_long %>%
  ggplot(data=.)+
  geom_raster(aes(x=Species,y=Proto_Karyotype,fill=Rank))+
  scale_fill_viridis(option="virdis") +
  geom_text(aes(x=Species,y=Proto_Karyotype,label=Rank),size=3,color="white")+
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust=0))+
  scale_x_discrete(position="top")+ylab("Proto Karyotype")

dev.off()
