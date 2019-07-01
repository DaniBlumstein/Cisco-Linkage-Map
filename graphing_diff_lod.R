#install.packages('reshape')
library(reshape)
#library(dplyr)
library(reshape2)
library(ggplot2)

new_data <- read.table("F:/LOD_parameter_space/consensus/consensus_hap_dip_lepmap_separate9.txt", quote="\"")
base_data <- read_excel("F:/LOD_parameter_space/consensus/consensus_hap_dip_separate_base.xlsx")

base_data$LG <- new_data$V1

base_data_t = dcast(base_data, LG ~ locus_type, value.var = "locus_type")
base_data_t

base_data_t_m = melt(base_data_t, id.vars = "LG" )

bar_graph <- ggplot(data =base_data_t_m, aes(LG, value, fill=variable)) +
  geom_bar(stat="identity", position = "stack") +
  ggtitle("LOD9 SL100")
bar_graph

