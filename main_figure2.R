
#pooja.singh09@gmail.com
#Singh et al L.callipterus 


library(rcompanion)
library(FSA)
library(ggplot2)
library(psych)


################################## low cov s23 < 2.4Mbp #############################


bm <- read.table("bm.allscaffs.snpden.txt", header=T)
dm <- read.table("dm.allscaffs.snpden.txt", header=T)
f <- read.table("females.allscaffs.snpden.txt", header=T)


dm_mean <- mean(dm$SNP_COUNT)
bm_mean <- mean(bm$SNP_COUNT)
f_mean <-mean(f$SNP_COUNT)



########### BM
### filling x ###


cov <- read.table("genome_wide_cov_MF_subset.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  unique(cov1$scaffold)
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input1 <- data.frame(matrix(NA, nrow = 50, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold == scaffs[i],]
	summ <- Summarize(dat$log2_bM.F)
	input1[i,1] <- paste("scaffold_", scaffs[i], sep="")
	input1[i,2] <- summ[5]
	input1[i,3] <- summ[7]
	input1[i,4] <- summ[2]
	input1[i,5] <- summ[3]
	
}

colnames(input1) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")




### filling y ###


bm <- read.table("bm.allscaffs.snpden_subset.txt", header=T)
dm <- read.table("dm.allscaffs.snpden_subset.txt", header=T)
f <- read.table("females.allscaffs.snpden", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))




input <- data.frame(matrix(NA, nrow = 50, ncol = 5))


for (i in 1:length(scaffs)){

	dat <- bmf[bmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/bm_mean
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/f_mean
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	input[i,1] <- scaffs[i]
	input[i,2] <- summ[5]
	input[i,3] <- summ[7]
	input[i,4] <- summ[2]
	input[i,5] <- summ[3]


}

colnames(input) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")

#### error bars bm cov and snp ####

input1 <-  input1[complete.cases(input1),]
input <-  input[complete.cases(input),]
input1 <- input1[rownames(input),]


colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey")
error.crosses(input1, input, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="bourgeois M:F log2 coverage", ylab="bourgeois M:F log2 SNP density", arrow.len=0)



################################## same cov s23 > 2.4Mbp #############################

############ BM
### filling x ###


cov <- read.table("genome_wide_cov_MF_subset1.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  unique(cov1$scaffold)
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input1a <- data.frame(matrix(NA, nrow = 50, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold2 == scaffs[i],]
	summ <- Summarize(dat$log2_bM.F)
	input1a[i,1] <- paste("scaffold_", scaffs[i], sep="")
	input1a[i,2] <- summ[5]
	input1a[i,3] <- summ[7]
	input1a[i,4] <- summ[2]
	input1a[i,5] <- summ[3]
	
}

colnames(input1a) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")




### filling y ###


bm <- read.table("bm.allscaffs.snpden_subset1.txt", header=T)
dm <- read.table("dm.allscaffs.snpden_subset1.txt", header=T)
f <- read.table("females.allscaffs.snpden", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))

inputa <- data.frame(matrix(NA, nrow = 50, ncol = 5))


for (i in 1:length(scaffs)){

	dat <- bmf[bmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/bm_mean
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/f_mean
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	inputa[i,1] <- scaffs[i]
	inputa[i,2] <- summ[5]
	inputa[i,3] <- summ[7]
	inputa[i,4] <- summ[2]
	inputa[i,5] <- summ[3]


}

colnames(inputa) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")

#### error bars bm cov and snp ####

inputa <-  inputa[complete.cases(inputa),]
input1a <-  input1a[complete.cases(input1a),]
input1a <- input1a[rownames(inputa),]

colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey")
error.crosses(input1a, inputa, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="bourgeois M:F log2 coverage", ylab="bourgeois M:F log2 SNP density", arrow.len=0)


################################## same cov s36 > 2 - 3.5 Mbp #############################

############ BM
### filling x ###


cov <- read.table("genome_wide_cov_MF_subset2.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  "36"
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input7a <- data.frame(matrix(NA, nrow = 1, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold == scaffs[i],]
	summ <- Summarize(dat$log2_bM.F)
	input7a[i,1] <- paste("scaffold_36")
	input7a[i,2] <- summ[5]
	input7a[i,3] <- summ[7]
	input7a[i,4] <- summ[2]
	input7a[i,5] <- summ[3]
	
}

colnames(input7a) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")




### filling y ###


bm <- read.table("bm.allscaffs.snpden_subset2.txt", header=T)
#dm <- read.table("dm.allscaffs.snpden_subset2.txt", header=T)
f <- read.table("females.allscaffs.snpden.txt", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))

input8a <- data.frame(matrix(NA, nrow = 1, ncol = 5))


for (i in 1:length(scaffs)){

	dat <- bmf[bmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/bm_mean
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/f_mean
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	input8a[i,1] <- scaffs[i]
	input8a[i,2] <- summ[5]
	input8a[i,3] <- summ[7]
	input8a[i,4] <- summ[2]
	input8a[i,5] <- summ[3]


}

colnames(input8a) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")



################################## same cov s36 1-2 Mb and 3.5 - 4.9 Mbp #############################

############ BM
### filling x ###


cov <- read.table("genome_wide_cov_MF_subset3.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  "36"
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input9a <- data.frame(matrix(NA, nrow = 1, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold == scaffs[i],]
	summ <- Summarize(dat$log2_bM.F)
	input9a[i,1] <- paste("scaffold_36")
	input9a[i,2] <- summ[5]
	input9a[i,3] <- summ[7]
	input9a[i,4] <- summ[2]
	input9a[i,5] <- summ[3]
	
}

colnames(input9a) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")


### filling y ###


bm <- read.table("bm.allscaffs.snpden_subset3.txt", header=T)
#dm <- read.table("dm.allscaffs.snpden_subset3.txt", header=T)
f <- read.table("females.allscaffs.snpden.txt", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))

input10a <- data.frame(matrix(NA, nrow = 1, ncol = 5))


for (i in 1:length(scaffs)){

	dat <- bmf[bmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/bm_mean
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/f_mean
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	input10a[i,1] <- scaffs[i]
	input10a[i,2] <- summ[5]
	input10a[i,3] <- summ[7]
	input10a[i,4] <- summ[2]
	input10a[i,5] <- summ[3]


}

colnames(input10a) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")


#####################################################################################DM


################################## low cov s23 < 2Mbp #############################

########## DM

cov <- read.table("genome_wide_cov_MF_subset.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  unique(cov1$scaffold)
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input1 <- data.frame(matrix(NA, nrow = 50, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold == scaffs[i],]
	summ <- Summarize(dat$log2_dM.F)
	input1[i,1] <- paste("scaffold_", scaffs[i], sep="")
	input1[i,2] <- summ[5]
	input1[i,3] <- summ[7]
	input1[i,4] <- summ[2]
	input1[i,5] <- summ[3]
	
}

colnames(input1) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")

bm <- read.table("bm.allscaffs.snpden_subset.txt", header=T)
dm <- read.table("dm.allscaffs.snpden_subset.txt", header=T)
f <- read.table("females.allscaffs.snpden.txt", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))
input2 <- data.frame(matrix(NA, nrow = 50, ncol = 5))


scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))

for (i in 1:length(scaffs)){

	dat <- dmf[dmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/dm_mean
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/f_mean
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	input2[i,1] <- scaffs[i]
	input2[i,2] <- summ[5]
	input2[i,3] <- summ[7]
	input2[i,4] <- summ[2]
	input2[i,5] <- summ[3]


}

colnames(input2) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")



#### error bars bm cov and snp ####


input1 <-  input1[complete.cases(input1),]
input2 <-  input2[complete.cases(input2),]
input1 <- input1[rownames(input2),]


colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey")
error.crosses(input1, input2, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="dwarf M:F log2 coverage", ylab="dwarf M:F log2 SNP density", arrow.len=0)



################################## same cov s23 > 2Mbp #############################


######## DM

cov <- read.table("genome_wide_cov_MF_subset1.txt", header=T)
cov1 <- cov[cov$scaffold <50,]


scaffs <-  unique(cov1$scaffold)
#cov1$CHROM <- paste("scaffold_", cov1$scaffold, sep="")


input1a <- data.frame(matrix(NA, nrow = 50, ncol = 5))

for (i in 1:length(scaffs)){

	dat <- cov1[cov1$scaffold == scaffs[i],]
	summ <- Summarize(dat$log2_dM.F)
	input1a[i,1] <- paste("scaffold_", scaffs[i], sep="")
	input1a[i,2] <- summ[5]
	input1a[i,3] <- summ[7]
	input1a[i,4] <- summ[2]
	input1a[i,5] <- summ[3]
	
}

colnames(input1a) <- c("scaffold", "cov_Q1", "cov_Q3", "mean", "sd")


bm <- read.table("bm.allscaffs.snpden_subset1.txt", header=T)
dm <- read.table("dm.allscaffs.snpden_subset1.txt", header=T)
f <- read.table("females.allscaffs.snpden", header=T)
scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))
bmf <- merge(bm,f, by=c("CHROM", "BIN_START"))
dmf <- merge(dm,f, by=c("CHROM", "BIN_START"))

input2a <- data.frame(matrix(NA, nrow = 50, ncol = 5))


scaffs <- unique(paste("scaffold_", cov1$scaffold, sep=""))

for (i in 1:length(scaffs)){

	dat <- dmf[dmf$CHROM == scaffs[i],]
	dat$SNP_COUNT.x.norm <- dat$SNP_COUNT.x/mean(dat$SNP_COUNT.x)
	dat$SNP_COUNT.y.norm <- dat$SNP_COUNT.y/mean(dat$SNP_COUNT.y)
	dat$logmf <- log(dat$SNP_COUNT.x.norm/dat$SNP_COUNT.y.norm)
	dat <- na.omit(dat)
	summ <- Summarize(dat$logmf)
	input2a[i,1] <- scaffs[i]
	input2a[i,2] <- summ[5]
	input2a[i,3] <- summ[7]
	input2a[i,4] <- summ[2]
	input2a[i,5] <- summ[3]


}

colnames(input2a) <- c("scaffold", "snpden_Q1", "snpden_Q3", "mean", "sd")


#### error bars bm cov and snp ####



input1a <-  input1a[complete.cases(input1a),]
input2a <-  input2a[complete.cases(input2a),]
input1a <- input1a[rownames(input2a),]



colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey")
error.crosses(input1a, input2a, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="dwarf M:F log2 coverage", ylab="dwarf M:F log2 SNP density", arrow.len=0)




######### plots dM

par <- input1a[input1a$scaffold == "scaffold_23",]
new_input_x <- rbind(input1, par)

par <- input2a[input2a$scaffold == "scaffold_23",]
new_input_y <- rbind(input2, par)

colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey", "grey","grey","grey","grey","grey", "grey","grey","grey","grey", "black")

error.crosses(new_input_x, new_input_y, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="dwarf M:F log2 coverage", ylab="dwarf M:F log2 SNP density", arrow.len=0, cex.lab=0.8, cex.axis=0.7)





######### plots bM


par <- inputa[inputa$scaffold == "scaffold_23",]
bm_new_input_y <- rbind(input, par)

par <- input8a[input8a$scaffold == "scaffold_36",]
bm_new_input_y1 <- rbind(bm_new_input_y, par)

par <- input10a[input10a$scaffold == "scaffold_36",]
bm_new_input_y2 <- rbind(bm_new_input_y1, par)


par <- input1a[input1a$scaffold == "scaffold_23",]
bm_new_input_x <- rbind(input1, par)

par <- input7a[input7a$scaffold == "scaffold_36",]
bm_new_input_x1 <- rbind(bm_new_input_x, par)

par <- input9a[input9a$scaffold == "scaffold_36",]
bm_new_input_x2 <- rbind(bm_new_input_x1, par)



bm_new_input_x2 <- bm_new_input_x2 [bm_new_input_x2$scaffold != "scaffold_7",]
bm_new_input_y2 <- bm_new_input_y2 [bm_new_input_y2$scaffold != "scaffold_4",]
bm_new_input_y2 <- bm_new_input_y2 [bm_new_input_y2$scaffold != "scaffold_20",]
bm_new_input_x2 <- bm_new_input_x2 [bm_new_input_x2$scaffold != "scaffold_32",]
dim(bm_new_input_y2)
dim(bm_new_input_x2)

bm_colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey", "grey","grey","grey","grey","grey","grey","grey","grey","black", "green", "blue")

error.crosses(bm_new_input_x2, bm_new_input_y2, sd=TRUE, labels=NA, colors=bm_colors, main=NA, xlab="bourgeois M:F log2 coverage", ylab="bourgeois M:F log2 SNP density", arrow.len=0, cex.lab=0.8, cex.axis=0.7)

### plot both


#pdf("main_figure_2_errorcross_2021.pdf",width=9, height=7)
svg("main_figure_2_errorcross_2021.svg",width=9, height=7)
par(mfrow=c(1,2))


bm_colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey", "grey","grey","grey","grey","grey","grey","grey","grey","black", "green", "blue")

error.crosses(bm_new_input_x2, bm_new_input_y2, sd=TRUE, labels=NA, colors=bm_colors, main=NA, xlab="bourgeois M:F log2 coverage", ylab="bourgeois M:F log2 SNP density", arrow.len=0, cex.lab=0.8, cex.axis=0.7, xlim=c(-0.3, 0.2), ylim=c(-0.3, 0.2))

abline(v=mean(bm_new_input_x2$mean), lty=2)
abline(h=mean(bm_new_input_y2$mean), lty=2)
legend(-0.3, -0.22, legend=c("scaffold_23  1 bp - 2.4 Mb", "scaffold_23  2.4 - 6.8 Mb", "scaffold_36  2.0 - 3.5 Mb", "scaffold_36  1.0 - 2.0 Mb & 3.5 - 4.9 Mb", "other scaffolds"),col=c("deeppink", "black","green", "blue", "grey"), lty=1, cex=0.7, box.lty=0, lwd=2)



colors <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","deeppink", "grey","grey","grey","grey","grey", "grey","grey","grey","grey","grey", "grey","grey","grey","grey", "black")

error.crosses(new_input_x, new_input_y, sd=TRUE, labels=NA, colors=colors, main=NA, xlab="dwarf M:F log2 coverage", ylab="dwarf M:F log2 SNP density", arrow.len=0, cex.lab=0.8, cex.axis=0.8, xlim=c(-0.3, 0.2), ylim=c(-0.3, 0.2))

abline(v=mean(new_input_x$mean), lty=2)
abline(h=mean(new_input_y$mean), lty=2)

legend(-0.3, -0.23, legend=c("scaffold_23  1 bp - 2.4 Mb", "scaffold_23  2.4 - 6.8 Mb", "other scaffolds"),col=c("deeppink", "black","grey"), lty=1, cex=0.7, box.lty=0, lwd=2)

dev.off()






