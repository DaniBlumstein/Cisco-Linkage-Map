rm(list = ls())
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
install.packages("qtl2convert", repos="http://rqtl.org/qtl2cran")
library(devtools)
install_github("kbroman/qtl")
library(qtl2, qtl2convert, qtl)


#bringing in the data 
P7_18 <- qtl::read.cross("csv", "F:/rqtl2/", "18_P7_rqtl_input.csv",  genotypes=c("AA","AB", "BB"))
jittered_P7_18 <- qtl::jittermap(P7_18)
P7_18_convereted <- qtl2::convert2cross2(P7_18)

#Calculating genotype probabilities
map <- qtl2::insert_pseudomarkers(P7_18_convereted$gmap, step=1)
pr <- qtl2::calc_genoprob(P7_18_convereted, map, error_prob=0.002)
apr <- qtl2::genoprob_to_alleleprob(pr)

#Calculating a kinship matrix
kinship <- qtl2::calc_kinship(pr, omit_x=TRUE)

#one may wish to eliminate the effect of varying marker density across the genome, and only use the probabilities along the grid of pseudomarkers (defined by the step argument to insert_pseudomarkers()). To do so, we need to first use calc_grid() to determine the grid of pseudomarkers, and then probs_to_grid() to omit probabilities for positions that are not on the grid.
grid <- qtl2::calc_grid(P7_18_convereted$gmap, step=1)
pr_grid <- qtl2::probs_to_grid(pr, grid)
kinship_grid <- qtl2::calc_kinship(pr_grid)

#Performing a genome scan by Haley-Knott regression
#The output of scan1() is a matrix of LOD scores, positions × phenotypes.
out <- qtl2::scan1(pr, P7_18_convereted$pheno)

#graphing LOD scores for each trait on each chr
par(mfrow = c(2, 2))
ymx <- qtl2::maxlod(out)
plot(out, map, lodcolumn=1, col="blue", ylim=c(0, ymx*1.02))
abline(h=3, col="red")
title("total length")

plot(out, map, lodcolumn=2, col="violetred", ylim=c(0, ymx*1.02))
title("standard length")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="purple", ylim=c(0, ymx*1.02))
title("body depth")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="orange", ylim=c(0, ymx*1.02))
title("weight")
abline(h=3, col="red")

title("18_P7", line = -1, outer = TRUE)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)


#used to derive the LOD support or Bayes credible intervals for QTL
qtl2::bayes_int(out, map, lodcolumn=2, chr=12, prob=0.95)

#Performing a genome scan with a linear model
out_pg <- qtl2::scan1(pr, P7_18_convereted$pheno, kinship)

#Here is a plot of the LOD scores, by Haley-Knott regression and the linear mixed model using either the standard kinship matrix or the LOCO method.
#standard length
color <- c("slateblue", "violetred")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(qtl2::maxlod(out), qtl2::maxlod(out_pg))
for(i in 1:4) {
  plot(out, map, lodcolumn=i, col=color[1], main=colnames(P7_18_convereted$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
  legend("topleft", lwd=2, col=color, c("H-K", "LMM"), bg="gray90", lty=c(1,1,2))
}

#Performing a permutation test
operm <- qtl2::scan1perm(pr, P7_18_convereted$pheno, n_perm=1000)
summary(operm)
summary(operm, alpha=c(0.2, 0.05))

#Esimated QTL effects
c2eff <- qtl2::scan1coef(pr[,"12"], P7_18_convereted$pheno[,"total_length"])

#plotting the effects
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["12"], columns=1:2, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

#Finally, plot the raw phenotypes against the genotypes at a single putative QTL position
g <- qtl2::maxmarg(pr, map, chr=20, pos=20.813, return_char=TRUE)

par(mar=c(4.1, 4.1, 0.6, 0.6))
qtl2::plot_pxg(g, P7_18_convereted$pheno[,"standard_length"], ylab="standard_length phenotype")


#################################################################################################
rm(list = ls())
#bringing in the data 
P7_24 <- qtl::read.cross("csv", "F:/rqtl2/", "24_P7_rqtl_input.csv",  genotypes=c("AA","AB", "BB"))
P7_24_convereted <- qtl2::convert2cross2(P7_24)

#Calculating genotype probabilities
map <- qtl2::insert_pseudomarkers(P7_24_convereted$gmap, step=1)
pr <- qtl2::calc_genoprob(P7_24_convereted, map, error_prob=0.002)
apr <- qtl2::genoprob_to_alleleprob(pr)

#Calculating a kinship matrix
kinship <- qtl2::calc_kinship(pr, omit_x=TRUE)

#one may wish to eliminate the effect of varying marker density across the genome, and only use the probabilities along the grid of pseudomarkers (defined by the step argument to insert_pseudomarkers()). To do so, we need to first use calc_grid() to determine the grid of pseudomarkers, and then probs_to_grid() to omit probabilities for positions that are not on the grid.
grid <- qtl2::calc_grid(P7_24_convereted$gmap, step=1)
pr_grid <- qtl2::probs_to_grid(pr, grid)
kinship_grid <- qtl2::calc_kinship(pr_grid)

#Performing a genome scan by Haley-Knott regression
#The output of scan1() is a matrix of LOD scores, positions × phenotypes.
out <- qtl2::scan1(pr, P7_24_convereted$pheno)

#graphing LOD scores for each trait on each chr
par(mfrow = c(2, 2))
ymx <- qtl2::maxlod(out)

plot(out, map, lodcolumn=1, col="blue", ylim=c(0, ymx*1.02))
abline(h=3, col="red")
title("total length")

plot(out, map, lodcolumn=2, col="violetred", ylim=c(0, ymx*1.02))
title("standard length")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="purple", ylim=c(0, ymx*1.02))
title("body depth")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="yellow", ylim=c(0, ymx*1.02))
title("weight")
abline(h=3, col="red")

title("24_P7", line = -1, outer = TRUE)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)

#####
# lodindex       lodcolumn chr    pos      lod  ci_lo  ci_hi
# 1        1    total_length  18 40.000 3.167314 32.275 47.206
# 2        1    total_length  29 32.958 3.395242 31.479 35.429
# 3        2 standard_length  29 32.958 3.468991 31.479 35.429
# 4        3      body_depth  18 40.000 3.648918 38.118 41.625
# 5        3      body_depth  29 32.958 3.109636 31.479 35.429

#used to derive the LOD support or Bayes credible intervals for QTL
qtl2::bayes_int(out, map, lodcolumn=2, chr=29, prob=0.95)

#Performing a genome scan with a linear mixed model
out_pg <- qtl2::scan1(pr, P7_24_convereted$pheno, kinship)

#Here is a plot of the LOD scores, by Haley-Knott regression and the linear mixed model using either the standard kinship matrix or the LOCO method.
#standard length
color <- c("slateblue", "violetred")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(qtl2::maxlod(out), qtl2::maxlod(out_pg))
for(i in 1:4) {
  plot(out, map, lodcolumn=i, col=color[1], main=colnames(P7_24_convereted$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
  legend("topleft", lwd=2, col=color, c("H-K", "LMM"), bg="gray90", lty=c(1,1,2))
}

#Performing a permutation test
operm <- qtl2::scan1perm(pr, P7_24_convereted$pheno, n_perm=1000)
summary(operm)
summary(operm, alpha=c(0.2, 0.05))

#Esimated QTL effects
c2eff <- qtl2::scan1coef(pr[,"29"], P7_24_convereted$pheno[,"total_length"])

#plotting the effects
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["12"], columns=1:2, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

#Finally, plot the raw phenotypes against the genotypes at a single putative QTL position
g <- qtl2::maxmarg(pr, map, chr=20, pos=20.813, return_char=TRUE)

par(mar=c(4.1, 4.1, 0.6, 0.6))
qtl2::plot_pxg(g, P7_24_convereted$pheno[,"standard_length"], ylab="standard_length phenotype")


#################################################################################################
rm(list = ls())
#bringing in the data 
P2_24 <- qtl::read.cross("csv", "F:/rqtl2/", "24_P2_rqtl_input.csv",  genotypes=c("AA","AB", "BB"))
P2_24_convereted <- qtl2::convert2cross2(P2_24)

#Calculating genotype probabilities
map <- qtl2::insert_pseudomarkers(P2_24_convereted$gmap, step=1)
pr <- qtl2::calc_genoprob(P2_24_convereted, map, error_prob=0.002)
apr <- qtl2::genoprob_to_alleleprob(pr)

#Calculating a kinship matrix
kinship <- qtl2::calc_kinship(pr, omit_x=TRUE)

#one may wish to eliminate the effect of varying marker density across the genome, and only use the probabilities along the grid of pseudomarkers (defined by the step argument to insert_pseudomarkers()). To do so, we need to first use calc_grid() to determine the grid of pseudomarkers, and then probs_to_grid() to omit probabilities for positions that are not on the grid.
grid <- qtl2::calc_grid(P2_24_convereted$gmap, step=1)
pr_grid <- qtl2::probs_to_grid(pr, grid)
kinship_grid <- qtl2::calc_kinship(pr_grid)

#Performing a genome scan by Haley-Knott regression
#The output of scan1() is a matrix of LOD scores, positions × phenotypes.
out <- qtl2::scan1(pr, P2_24_convereted$pheno)

#graphing LOD scores for each trait on each chr
par(mfrow = c(2, 2))
ymx <- qtl2::maxlod(out)

plot(out, map, lodcolumn=1, col="blue", ylim=c(0, ymx*1.02))
abline(h=3, col="red")
title("total length")

plot(out, map, lodcolumn=2, col="violetred", ylim=c(0, ymx*1.02))
title("standard length")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="purple", ylim=c(0, ymx*1.02))
title("body depth")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="yellow", ylim=c(0, ymx*1.02))
title("weight")
abline(h=3, col="red")

title("24_P2", line = -1, outer = TRUE)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)

# lodindex  lodcolumn chr    pos      lod ci_lo  ci_hi
# 1        3 body_depth   6 31.786 4.805940 2.036 42.409
# 2        4     weight   6 31.786 3.867253 2.036 42.409

#used to derive the LOD support or Bayes credible intervals for QTL
qtl2::bayes_int(out, map, lodcolumn=2, chr=12, prob=0.95)

#Performing a genome scan with a linear mixed model
out_pg <- qtl2::scan1(pr, P2_24_convereted$pheno, kinship)

#Here is a plot of the LOD scores, by Haley-Knott regression and the linear mixed model using either the standard kinship matrix or the LOCO method.
#standard length
color <- c("slateblue", "violetred")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(qtl2::maxlod(out), qtl2::maxlod(out_pg))
for(i in 1:4) {
  plot(out, map, lodcolumn=i, col=color[1], main=colnames(P2_24_convereted$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
  legend("topleft", lwd=2, col=color, c("H-K", "LMM"), bg="gray90", lty=c(1,1,2))
}

#Performing a permutation test
operm <- qtl2::scan1perm(pr, P2_24_convereted$pheno, n_perm=1000)
summary(operm)
summary(operm, alpha=c(0.2, 0.05))

#Esimated QTL effects
c2eff <- qtl2::scan1coef(pr[,"12"], P2_24_convereted$pheno[,"total_length"])

#plotting the effects
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["12"], columns=1:2, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

#Finally, plot the raw phenotypes against the genotypes at a single putative QTL position
g <- qtl2::maxmarg(pr, map, chr=20, pos=20.813, return_char=TRUE)

par(mar=c(4.1, 4.1, 0.6, 0.6))
qtl2::plot_pxg(g, P2_24_convereted$pheno[,"standard_length"], ylab="standard_length phenotype")


#################################################################################################
rm(list = ls())
#bringing in the data 
P9_24 <- qtl::read.cross("csv", "F:/rqtl2/", "24_P9_rqtl_input.csv",  genotypes=c("AA","AB", "BB"))
P9_24_convereted <- qtl2::convert2cross2(P9_24)

#Calculating genotype probabilities
map <- qtl2::insert_pseudomarkers(P9_24_convereted$gmap, step=1)
pr <- qtl2::calc_genoprob(P9_24_convereted, map, error_prob=0.002)
apr <- qtl2::genoprob_to_alleleprob(pr)

#Calculating a kinship matrix
kinship <- qtl2::calc_kinship(pr, omit_x=TRUE)

#one may wish to eliminate the effect of varying marker density across the genome, and only use the probabilities along the grid of pseudomarkers (defined by the step argument to insert_pseudomarkers()). To do so, we need to first use calc_grid() to determine the grid of pseudomarkers, and then probs_to_grid() to omit probabilities for positions that are not on the grid.
grid <- qtl2::calc_grid(P9_24_convereted$gmap, step=1)
pr_grid <- qtl2::probs_to_grid(pr, grid)
kinship_grid <- qtl2::calc_kinship(pr_grid)

#Performing a genome scan by Haley-Knott regression
#The output of scan1() is a matrix of LOD scores, positions × phenotypes.
out <- qtl2::scan1(pr, P9_24_convereted$pheno)

#graphing LOD scores for each trait on each chr
par(mfrow = c(2, 2))
ymx <- qtl2::maxlod(out)

plot(out, map, lodcolumn=1, col="blue", ylim=c(0, ymx*1.02))
abline(h=3, col="red")
title("total length")

plot(out, map, lodcolumn=2, col="violetred", ylim=c(0, ymx*1.02))
title("standard length")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="purple", ylim=c(0, ymx*1.02))
title("body depth")
abline(h=3, col="red")

plot(out, map, lodcolumn=2, col="yellow", ylim=c(0, ymx*1.02))
title("weight")
abline(h=3, col="red")

title("24_P9", line = -1, outer = TRUE)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)
#finding LOD peaks
qtl2::find_peaks(out, map, threshold=3, drop=1.5)

#used to derive the LOD support or Bayes credible intervals for QTL
qtl2::bayes_int(out, map, lodcolumn=2, chr=12, prob=0.95)

#Performing a genome scan with a linear mixed model
out_pg <- qtl2::scan1(pr, P9_24_convereted$pheno, kinship)

#Here is a plot of the LOD scores, by Haley-Knott regression and the linear mixed model using either the standard kinship matrix or the LOCO method.
#standard length
color <- c("slateblue", "violetred")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(qtl2::maxlod(out), qtl2::maxlod(out_pg))
for(i in 1:4) {
  plot(out, map, lodcolumn=i, col=color[1], main=colnames(P9_24_convereted$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
  legend("topleft", lwd=2, col=color, c("H-K", "LMM"), bg="gray90", lty=c(1,1,2))
}

#Performing a permutation test
operm <- qtl2::scan1perm(pr, P9_24_convereted$pheno, n_perm=1000)
summary(operm)
summary(operm, alpha=c(0.2, 0.05))

#Esimated QTL effects
c2eff <- qtl2::scan1coef(pr[,"12"], P9_24_convereted$pheno[,"total_length"])

#plotting the effects
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["12"], columns=1:2, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

#Finally, plot the raw phenotypes against the genotypes at a single putative QTL position
g <- qtl2::maxmarg(pr, map, chr=20, pos=20.813, return_char=TRUE)

par(mar=c(4.1, 4.1, 0.6, 0.6))
qtl2::plot_pxg(g, P9_24_convereted$pheno[,"standard_length"], ylab="standard_length phenotype")
