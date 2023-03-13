

#pooja.singh09@gmail.com
#Singh et al Dwarfism L. calllipterus
#Figure 1A - D


library(scales)
library(Rmisc)

#### genome wide cov plot for males versus females

a <- read.table("cov_norm_genomewide_141scaffs_2021.txt", sep='\t',header=T)

dm <- CI(a$log2_dM.F, ci=0.99)
dm_mean= dm[[2]]
dm_upperCI= dm[[1]]
dm_lowerCI= dm[[3]]
bm <- CI(a$log2_bM.F, ci=0.99)
bm_mean= bm[[2]]
bm_upperCI= bm[[1]]
bm_lowerCI= bm[[3]]

dq <- quantile(a$log2_dM.F, probs=c(0.01, 0.99))
bq <- quantile(a$log2_bM.F, probs=c(0.01, 0.99))

########################################## main fig 1

#pdf("main_figure_1.pdf",width=15, height=7)
svg("main_figure_1_revision2021.svg",width=15, height=7)
line = 1
cex = 1.5
side = 3
adj=-0.05

palette(c("grey","black"))

tick=c(1,49,89,129,163,195,226,257,285,313,338,361,384,406,427,448,468,487,505,521,537,552,566,580,594,608,622,636,649,662,675,688,700,712,724,735,746,756,766,776,786,796,806,816,826,836,846,856,866,875,884,893,902,911,920,929,938,947,956,965,973,981,989,997,1005,1013,1021,1029,1037,1045,1052,1059,1066,1073,1080,1087,1094,1101,1108,1114,1120,1126,1132,1138,1144,1150,1156,1162,1168,1174,1180,1186,1192,1198,1203,1208,1213,1218,1223,1228,1233,1238,1243,1248,1253,1258,1263,1268,1273,1278,1283,1288,1292,1296,1300,1304,1308,1312,1316,1320,1324,1328,1332,1336,1340,1344,1348,1352,1356,1360,1364,1368,1372,1376,1380,1384,1388,1392,1396,1400,1404,1408)
#tick1 = c(25,69,109,146,179,210.5,241.5,271,299,325.5,349.5,372.5,395,416.5,437.5,459,478.5,496,513,529,544.5,559,573,587,608,629,642.5,655.5,668.5,681,700)
#par(mfrow=c(2,1), oma = c(1,4,2,0) + 0.1,mar = c(2,4,0,0) + 1.5)

par(mar=c(4,4,3,2))
layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(2, 2), # Heights of the two rows
       widths = c(2, 1))



plot(a$number, a$log2_bM.F , col="grey", pch=20, cex=0.5, ylim=c(-0.25, 0.3), xaxt='n', ylab="bourgeois M:F log2 coverage", cex.lab=1,font=1,xlab="scaffold")
lines(a$avrg_bM, col="blue", lwd=2)
abline(h=bm_mean, lty=2)
#abline(h=bm_upperCI, lty=2, col="grey")
#abline(h=bm_lowerCI, lty=2, col="grey")
abline(h=bq[1], lty=2, col="grey")
abline(h=bq[2], lty=2, col="grey")

axis(side = 1, at = tick, labels=NA, cex.axis = 0.8)
#axis(side = 1, at = tick1, tick=F, labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), cex.axis = 0.5,font=2)
rect(577,-0.3,594,0.32, col= rgb(0,0,1.0,alpha=0.2), border=NA)
#mtext(side=1, line=3, "scaffold", cex=1,font=1)
mtext("A", side=side, line=line, cex=cex, adj=adj,font=2)

s23 <- a[a$scaffold == "23",]
plot(density(s23$log2_bM.F), col="deeppink", xlim=c(-0.32,0.1), ylim=c(0,55), main=NA, xlab="bourgeois M:F log2 coverage", ylab="Density")
legend(-0.3, 50, legend=c("scaffold_23", "scaffold_36", "other scaffolds"),col=c("deeppink", "green", "grey"), lty=1, cex=1, box.lty=0, lwd=2)
mtext("B", side=side, line=line, cex=cex, adj=adj,font=2)
s36 <- a[a$scaffold == "36",]
lines(density(s36$log2_bM.F),col="green")

scaffs<- seq(from = 0, to = 40)
scaffs<- scaffs[-24];
scaffs<- scaffs[-37];
for (i in scaffs){

	dat <- a[a$scaffold == i,]
	lines(density(dat$log2_bM.F),col=alpha("grey", 0.4))
}



plot(a$number, a$log2_dM.F, col="grey", pch=20, cex=0.5, ylim=c(-0.25, 0.3), xaxt='n', ylab="dwarf M:F log2 coverage", cex.lab=1,font=1,xlab="scaffold")
lines(a$avrg_dM, col="red", lwd=2)
abline(h=dm_mean, lty=2)
#abline(h=dm_upperCI, lty=2, col="grey")
#abline(h=dm_lowerCI, lty=2, col="grey")
abline(h=dq[1], lty=2, col="grey")
abline(h=dq[2], lty=2, col="grey")
axis(side = 1, at = tick, labels=NA, cex.axis = 0.8)
#axis(side = 1, at = tick1, tick=F, labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30), cex.axis = 0.5,font=2)
rect(577,-0.3,594,0.32, col= rgb(1,0,0,alpha=0.2), border=NA)
#mtext(side=1, line=3, "scaffold", cex=1,font=1)
mtext("C", side=side, line=line, cex=cex, adj=adj,font=2)

plot(density(s23$log2_dM.F), col="deeppink", xlim=c(-0.32,0.1), ylim=c(0,55), main=NA, xlab="dwarf M:F log2 coverage", ylab="Density")
legend(-0.3, 50, legend=c("scaffold_23", "other scaffolds"),col=c("deeppink", "grey"), lty=1, cex=1, box.lty=0, lwd=2)
mtext("D", side=side, line=line, cex=cex, adj=adj,font=2)
scaffs<- seq(from = 0, to = 40)
scaffs<- scaffs[-24];
for (i in scaffs){

	dat <- a[a$scaffold == i,]
	lines(density(dat$log2_dM.F),col=alpha("grey", 0.4))
}

dev.off()
