#pooja.singh09@gmail.com
#Singh et al L.callipterus 

library(ggplot2)
library(zoo)
library(cowplot)
require(gridExtra)
library(Rmisc)

#### genome wide cov plot for males versus females

a <- read.table("genome_wide_cov_MF.txt", sep='\t',header=T)

dm <- CI(a$log2_dM.F, ci=0.95)
dm_mean= dm[[2]]
dm_upperCI= dm[[1]]
dm_lowerCI= dm[[3]]
bm <- CI(a$log2_bM.F, ci=0.95)
bm_mean= bm[[2]]
bm_upperCI= bm[[1]]
bm_lowerCI= bm[[3]]

dq <- quantile(a$log2_dM.F, probs=c(0.05, 0.95))
bq <- quantile(a$log2_bM.F, probs=c(0.05, 0.95))

######



#### s23 coverage #####


c <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020/s23_cov.txt", header=T)

plot1  <- ggplot(c) + ylim(-0.6, 0.1) + geom_point(aes(x=pos_start, y=log_bM.F), col="grey") + geom_line(aes(x=pos_start, y=log_bM.F), col="blue") + theme_bw() + xlab("scaffold_23 position (bp)") + ylab("bourgeois M:F \n log2 coverage") + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=bq[1], linetype="dashed", color = "grey") + geom_hline(yintercept=bq[2], linetype="dashed", color = "grey") + geom_hline(yintercept=bm_mean, linetype="dashed")


plot2 <- ggplot(c)  + ylim(-0.6, 0.1)  + geom_point(aes(x=pos_start, y=log_dM.F), col="grey") + geom_line(aes(x=pos_start, y=log_dM.F), col="red") + theme_bw() + xlab("scaffold_23 position (bp)") + ylab("dwarf M:F \n log2 coverage")  + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=dq[1], linetype="dashed", color = "grey") + geom_hline(yintercept=dq[2], linetype="dashed", color = "grey") + geom_hline(yintercept=dm_mean, linetype="dashed")

## snp density ######


bm <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/snpdensity/all.filt.vcf.s23.bm.SNP.recode.vcf.30kb.snpden", header=T)
dm <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/snpdensity/all.filt.vcf.s23.dm.SNP.recode.vcf.30kb.snpden", header=T)
f <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/snpdensity/all.filt.vcf.s23.females.SNP.recode.vcf.30kb.snpden", header=T)

dm_mean <- mean(dm$SNP_COUNT)
bm_mean <- mean(bm$SNP_COUNT)
f_mean <-mean(f$SNP_COUNT)

dmf <- merge(dm,f, by="BIN_START")
bmf <- merge(bm,f, by="BIN_START")
dmf$mf <- (dmf$SNP_COUNT.x/dm_mean)/(dmf$SNP_COUNT.y/f_mean)
bmf$mf <- (bmf$SNP_COUNT.x/bm_mean)/(bmf$SNP_COUNT.y/f_mean)


plot3 <- ggplot(bmf) + geom_point(aes(x=BIN_START, y=log(mf)), col="grey") + geom_line(aes(x=BIN_START, y=log(rollmean(mf, 10, na.pad=TRUE))), col="blue") + theme_bw() + xlab("scaffold_23 position (bp)") + ylab("bourgeois M:FSNP density")  + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + ylim(-0.5,0.2) + geom_hline(yintercept=-0.12501439, linetype="dashed",  col="grey") + geom_hline(yintercept=0.09733678, linetype="dashed", col="grey") + geom_hline(yintercept=bm_mean, linetype="dashed") + geom_hline(yintercept=0.0009056152, linetype="dashed") 

plot4 <- ggplot(dmf) + geom_point(aes(x=BIN_START, y=log(mf)), col="grey") + geom_line(aes(x=BIN_START, y=log(rollmean(mf, 10, na.pad=TRUE))), col="red") + theme_bw() + xlab("scaffold_23 position (bp)") + ylab("dwarf M:F SNP density")  + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + ylim(-0.5,0.2) + geom_hline(yintercept=-0.2050186, linetype="dashed",  col="grey") + geom_hline(yintercept=0.1547420, linetype="dashed", col="grey") + geom_hline(yintercept=dm_mean, linetype="dashed")+ geom_hline(yintercept=-0.002519406, linetype="dashed") 


#### s23 FST males versus femles################   

all <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/popgenome/all_popgenome_30000.txt", header=T)
s23 <- all[all$chr == "scaffold_23", ]
bf <- s23[s23$pop == "bm",]
df <- s23[s23$pop == "dm",]
dm_mean <- mean(df$FST)
bm_mean <- mean(bf$FST)

stats <- read.table("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/popgenome/all_popgenome_30000_CI_quantiles.txt", header=T)

plot5 <-ggplot(bf) + geom_point(aes(x=start, y=FST), col="grey")+ ylim(0,0.5)  + theme_bw()+ geom_line(aes(x=start, y=rollmean(FST,10, na.pad=TRUE)), col="blue") + xlab("scaffold_23 position (bp)") + ylab("bourgeois M vs F FST")  + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=quantile(bf$FST, probs=c(0.05))[[1]], linetype="dashed", col="grey")  + geom_hline(yintercept=quantile(bf$FST, probs=c(0.95))[[1]], linetype="dashed", col="grey") + geom_hline(yintercept=bm_mean, linetype="dashed")

plot6 <-ggplot(df) + geom_point(aes(x=start, y=FST), col="grey")+ ylim(0,0.5) + theme_bw() + geom_line(aes(x=start, y=rollmean(FST,10, na.pad=TRUE)), col="red") + xlab("scaffold_23 position (bp)") + ylab("dwarf M vs F FST")  + theme(axis.title=element_text(size=10)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=quantile(df$FST, probs=c(0.05))[[1]], linetype="dashed", col="grey")  + geom_hline(yintercept=quantile(df$FST, probs=c(0.95))[[1]], linetype="dashed", col="grey") + geom_hline(yintercept=dm_mean, linetype="dashed") 



pdf("main_figure_3.pdf",width=9, height=10)
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2,nrow=3, align="v")
dev.off()

svg("main_figure_3.svg",width=9, height=10)
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2,nrow=3, align="v")
dev.off()
