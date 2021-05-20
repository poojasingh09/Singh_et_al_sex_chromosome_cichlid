# pooja.singh09@gmail.com
# jan20201
# Script for running PopGen data for L. callipterus sex chr Paper - Singh et al.
# Calculates male/female Fst, dxy and pi for each population
# Caldulates Da sensu Fraser et al 2020 (M-F Dxy - F Pi)
# Calculates these for 30kb windows
##################################################################

# Set
setwd("/Volumes/pooja2019/PhD_2015/PhD_research_projects/Dwarfism/2020sept/snps/popgenome")

# Get libs
lib<-as.vector(c("pbapply","PopGenome","ggplotify","cowplot","gridExtra","viridis","grid","data.table","HiTC","Sushi","ggplot2","parallel","dplyr","tidyr"))
lapply(lib,library,character.only=TRUE)

vcf<-"all.filt.vcf.s23.callipterus.SNP.recode.vcf.gz"

# Catch errors
wind_size=30000


# Get length
chr_length <- read.table("N_brichardi_v1.assembly.fasta.fai.s23", header=F)
chrs<- data.frame(chr_length[,1])
tid <- chrs[1,]
frompos <- 1
topos <- chr_length$V2

# Read in the VCF
snp <- readVCF(file = vcf, tid=tid, frompos = 1, topos = topos, numcols=1000000, include.unknown=TRUE)

# Set the males and females
c_M = c("bourgeois_male_1","bourgeois_male_2","dwarf_male_1","dwarf_male_2")
c_dM = c("dwarf_male_1","dwarf_male_2")
c_bM = c("bourgeois_male_1","bourgeois_male_2")
c_F = c("bourgeois_male_1_daughter","bourgeois_male_2_daughter","dwarf_male_1_daughter","dwarf_male_2_daughter")

c_list<-list(c_M,c_F)
d_list<-list(c_dM,c_F)
b_list<-list(c_bM,c_F)


comparison_lists<-list(c_list, b_list, d_list, m_list)
names<-c("cal", "bm", "dm")
 
# Lapply over list of comparisons
chr_popgen<-data.frame(rbindlist(mclapply(1:length(comparison_lists),function(x){

snp<-set.populations(snp,do.call("list",comparison_lists[[x]]),diploid=TRUE)
snp<-set.outgroup(snp,FALSE)

# Split into windows
win_SNP_30k<-sliding.window.transform(snp,width=wind_size,jump=wind_size,type=2)

#do pop stats
win_SNP_30k <-F_ST.stats(win_SNP_30k,mode="nucleotide")

win_SNP_30k <-neutrality.stats(win_SNP_30k,FAST=FALSE)

win_SNP_30k <-diversity.stats.between(win_SNP_30k,nucleotide.mode=TRUE)

# Get centre of the window
 genome.pos_30k <- sapply(win_SNP_30k@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
   val   <- mean(as.numeric(split))
    return(val)
  })
   
# Output results matrix
PG_out<-data.frame(chr = tid,
                   FST=win_SNP_30k@nucleotide.F_ST,
                   n.seg.M=win_SNP_30k@n.segregating.sites[,1],
                   n.seg.F=win_SNP_30k@n.segregating.sites[,2],
                   dxy=(win_SNP_30k@nuc.diversity.between/wind_size), 
                   pi.M=(win_SNP_30k@nuc.diversity.within/wind_size)[,1],
                   pi.F=(win_SNP_30k@nuc.diversity.within/wind_size)[,2],
                   start=genome.pos_30k -(wind_size/2)-1.5,
                   end=genome.pos_30k +(wind_size/2)-1.5,
                   pop=names[x])
colnames(PG_out)[c(2,5)]<-c("FST","DXY")
#PG_out<-na.omit(PG_out)
PG_out[PG_out$FST < 0,"FST"]<-0
rownames(PG_out)<-NULL

return(PG_out)
},mc.cores=3)))


## Calculate Da

chr_popgen$Da <- chr_popgen$DXY - chr_popgen$pi.F


write.table(na.omit(chr_popgen),
            paste0("all.filt.vcf.s23.callipterus.SNP.recode.vcf_",wind_size,"_",tid,".txt"),row.names=F,quote=F,sep="\t")


#plots

both <- chr_popgen[chr_popgen$pop == "cal",]
bm <- chr_popgen[chr_popgen$pop == "bm",]
dm <- chr_popgen[chr_popgen$pop == "dm",]


par(mfrow=c(3,1))
plot(both$start, both$DXY)
plot(both$start, both$Da)
plot(both$start, both$FST)


pdf("s23_dxy_da_fst.pdf")
par(mfrow=c(3,2))
plot(bm$start, bm$FST, ylim=c(0,0.5), xlab="position (bp)", pch=16, col="blue", ylab= "bourgeois males vs females FST")
abline(h=mean(bm$FST),lty=2)
plot(dm$start, dm$FST, ylim=c(0,0.5), xlab="position (bp)", pch=16, col="red", , ylab= "dwarf males vs females FST")
abline(h=mean(dm$FST),lty=2)
plot(bm$start, bm$DXY, ylim=c(0,0.006), xlab="position (bp)", pch=16, col="blue", ylab= "bourgeois males vs females DXY")
abline(h=mean(bm$DXY),lty=2)
plot(dm$start, dm$DXY, ylim=c(0,0.006), xlab="position (bp)", pch=16, col="red", ylab= "dwarf males vs females DXY")
abline(h=mean(dm$DXY),lty=2)
plot(bm$start, bm$Da, ylim=c(0,0.0017), xlab="position (bp)", pch=16, col="blue",ylab= "bourgeois males vs females Da")
abline(h=mean(bm$Da),lty=2)
plot(dm$start, dm$Da, ylim=c(0,0.0017), xlab="position (bp)", pch=16, col="red", ylab= "dwarf males vs females Da")
abline(h=mean(dm$Da),lty=2)
dev.off()



par(mfrow=c(3,1))
plot(both$start, both$pi.F)
plot(bm$start, bm$pi.M)
plot(dm$start, dm$pi.M)



pdf("s23_dxy_da_fst_2.pdf")
par(mfrow=c(3,1))
plot(both$start, both$FST, xlab="position (bp)", pch=16, col="black", ylab= "males vs females FST")
abline(h=mean(both$FST),lty=2)
plot(both$start, both$DXY, xlab="position (bp)", pch=16, col="grey", ylab= "males vs females DXY")
abline(h=mean(both$DXY),lty=2)
plot(both$start, both$Da, xlab="position (bp)", pch=16, col="pink", ylab= "males vs females Da")
abline(h=mean(both$Da),lty=2)

dev.off()





























