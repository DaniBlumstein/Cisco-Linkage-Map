library(circlize,RColorBrewer)
jpeg("circle_plot_color_v2.jpg",height=8,width=8,units="in",res=700)
#Customize plot
par(fig=c(0,1,0,1),mar=c(1,1,1,1))
LinkageMap<- read.delim("F:/LepMap3/consensus_files_for_lepmap/final_lepmap/consensus_hap_dip_lepmap3_lod15_final_renamed.txt")

circos.clear()

circos.par("track.height" = 0.1 , start.degree=95, cell.padding = c(0.005, 0, 0.005, 0), track.margin = c(0.0045, 0.0045))

circos.initialize(factors = LinkageMap$LG_NEW, x = LinkageMap$female_position)

circos.track(factors = LinkageMap$LG_NEW, y = LinkageMap$female_position,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, + uy(3, "mm"),
                           CELL_META$sector.index)
             })

#circos.trackPoints(LinkageMap$LG, LinkageMap$female_position, LinkageMap$male_position, col = LinkageMap$DuplicationStatus, pch = 16, cex = 0.5)
circos.link("1", c(36.87, 101.353), "38", c(0,58.97), col = "grey82",h = 0.65,border = "grey")

circos.link("3", c(0, 57.7), "27", c(0,57.996), col = "grey82",h = 0.65, border = "grey")

circos.link("6", c(40.8, 91.364), "18", c(0,62.155), col = "grey82",h = 0.65, border = "grey")

circos.link("7", c(44.88, 57.72), "28", c(0,57.793), col = "grey82",h = 0.65, border = "grey")

circos.link("16", c(0, 62.253), "33", c(0,56.371), col = "grey82",h = 0.65,border = "grey")

circos.link("17", c(0, 62.247), "32", c(0,56.557), col = "grey82",h = 0.65,border = "grey")


#AORes
#display.brewer.pal(n = 9, name = 'Blues')
brewer.pal(n = 8, name = "Blues")

circos.link("3", c(56.439, 92.033), "34", c(0.914,54.932), col = "#2171B5",h = 0.65)

circos.link("1", c(0, 37.235), "12", c(5.674,61.652), col = "#4292C6", h = 0.65)

circos.link("7", c(0, 43.952), "2", c(0.914,44.739), col = "#6BAED6", h = 0.65)

circos.link("8", c(0, 43.713), "2", c(51.252,96.059), col = "#9ECAE1", h = 0.65)

circos.link("5", c(0, 68.004), "9", c(26.571,63.916), col = "#C6DBEF", h = 0.65)

circos.link("6", c(0, 27.457), "15", c(4.361,24.539), col = "#DEEBF7", h = 0.65)

#higher support LORes
#display.brewer.pal(n = 9, name = 'Oranges')
brewer.pal(n = 9, name = "Oranges")

circos.link("29", c(25.627, 52.491), "37", c(46.031,49.663), col = "#A63603",h = 0.65)

circos.link("22", c(3.364, 12.893), "24", c(3.592,20.699), col = "#D94801", h = 0.65)

circos.link("4", c(52.996, 92.451), "8", c(47.341,55.252), col = "#F16913", h = 0.65)

circos.link("23", c(23.494, 56.153), "14", c(20.229,48.024), col = "#FD8D3C", h = 0.65)


#low support LORes
circos.link("21", c(2.591, 8.617), "30", c(3.934,8.085), col = "#FDAE6B",h = 0.65)

circos.link("25", c(3.573, 8.933), "26", c(1.315,6.827), col = "#FDAE6B",h = 0.65)

circos.link("5", c(80.456, 89.427), "36", c(44.669,48.669), col = "#FDAE6B",h = 0.65)

circos.link("10", c(63.491,67.491), "31", c(48.917,52.917), col = "#FDD0A2",h = 0.65)

circos.link("13", c(48.94,52.94), "11", c(53.582,57.582), col = "#FDD0A2",h = 0.65)

circos.link("19", c(13.86,17.86), "20", c(42.937,46.937), col = "#FDD0A2",h = 0.65)

circos.link("35", c(34.747,38.747), "4", c(5.367,9.367), col = "#FDD0A2",h = 0.65)




legend("bottomleft", inset=.005, bty = "n",
       c("AORe","LORe", "Theoretical"), fill=c("coral", rgb(0/255,112/255,192/255), "grey82"), border = c("black", "black", "black"), horiz=F, cex=0.8)

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'LORe')
text(x=1.5, y = seq(0,1,l=5), labels = c(42,46,52,58,69,86))
rasterImage(legend_image, 0, 0, 1,1)


dev.off()
# 
# lgd_lines = Legend(at = c("AORs", "LORs"), type = "lines", 
#                    legend_gp = gpar(col = 4:5, lwd = 2), title_position = "topleft", 
#                    title = "Track2")
# 
# lgd_list_vertical = packLegend(lgd_lines, direction = "horizontal")
# 
# par(omi = gridOMI(), new = TRUE)
# 
# draw(lgd_list_vertical, x = circle_size, just = "top")

# for (row in 1:(nrow(LinkageMap))){
#   e <- LinkageMap$LG[row]
#   d <- LinkageMap$female_position[row]
#   a <- LinkageMap$LG_friend[row]
#   b <- LinkageMap$Pos_friend[row]
#   fact <- as.factor(LinkageMap$LG_friend)[row]
#   circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.7, col = colvec[fact])
# }

