
# I looked at copy number and not named mutations (e.g. p.V216E for TP53)
# In the squam data I'm assuming blanks mean no amplification or deletion

# Hawazin noted:

# adeno have:
#   NKX2-1 amp
#   More frequent EGFR, KRAS mutations 
#   More frequent ALK fusions

# squam have:
#   SOX2, TP63, PIK3CA amp
#   -if no amp for these then NOTCH1, NOTCH2, or ASCL4  truncating mutations
#   More frequent TP53 mutations

# In the TCGA adeno paper download, I looked at "Focal Copy Number Alterations (N=no focal amplification or deletion; A=focal amplification; D=focal deletion)" for n=230 and ngenes=30.

# In the TCGA squam paper download, I looked at "Focal Copy Number Alterations (N=no focal amplification or deletion; A=focal amplification; D=focal deletion)" for n=178 and ngenes=43.

mydir<-"~/Documents/data/chucks5/D/3class/mutations"
setwd(mydir)

load("~/Documents/data/chucks5/D/3class/mutations/paper_downloads.RData")

#load("~/Documents/data/chucks5/D/3class/lsp_and_subtype/lsp_and_subtype.RData")
load("~/Documents/data/chucks5/D/3class/with_scc_lsp_and_subtype/lsp_and_subtype.RData")

### compare my scc subtype ###
tmp<-merge(sq.clin,scc.sub.tcga$myres,by=0)
addmargins(table(tmp$subtype,tmp$subtype.tcga,exclude=NULL)) # all 178 

### compare my ad subtype ###
tmp<-merge(ad.clin,sub.tcga$myres,by=0)
addmargins(table(tmp$subtype,tmp$expression_subtype,exclude=NULL)) # all 230

### compare clinical ###
common<-names(ad.clin)[names(ad.clin) %in% names(sq.clin)]
tmp<-merge(w.tcga,rbind(ad.clin[,common],sq.clin[,common]),by=0)
#plot(tmp$age,tmp$Age); abline(0,1)
table(tmp$T.stage.x,tmp$T.stage.y,exclude=NULL)
table(tmp$N.stage.x,tmp$N.stage.y,exclude=NULL)
table(tmp$sex,tmp$Sex,exclude=NULL)

###################
### copy number ###
###################

# 4 groups, plot proportion with Amp, Del, use all genes.  
# Indicate which genes were not in ad or not in sq data.
# -in squam only, or in adeno only (color the gene name).

ad.cn<-data.frame(ad.cn)
sq.cn<-data.frame(sq.cn)
ad.genes<-names(ad.cn)
sq.genes<-names(sq.cn)
genes<-sort(union(names(ad.cn),names(sq.cn)))

ad.cn[,genes[!genes %in% ad.genes]]<-NA
ad.cn<-ad.cn[,genes]
sq.cn[,genes[!genes %in% sq.genes]]<-NA
sq.cn<-sq.cn[,genes]

table(names(ad.cn)==names(sq.cn))
a<-rbind(ad.cn,sq.cn)
a<-data.frame(merge(a,lsp.tcga$l,by=0),row.names=1)

### table ###

sink("samplesize.txt") # note same for copy number and mutation data b/c all in one spreadsheet
addmargins(table(a$group,a$group.alt))
sink()

### p-values ###

# ad-notad vs other group (or ad-sq vs other group)?  2x2 table for Amp/not x g1/g2 and table for Del/not x g1/g2
#p.amp<-data.frame(gene=genes,vs.adad=NA,vs.sqsq=NA)
#p.del<-data.frame(gene=genes,vs.adad=NA,vs.sqsq=NA)

### proportions ###

prop.amp<-prop.del<-data.frame(gene=genes,adad=NA,adnotad=NA,sqsq=NA,sqnotsq=NA,row.names=1)
for(i in ad.genes){
  prop.amp[i,"adad"]<-mean(a[a$group=="ad-ad",i]=="A")
  prop.del[i,"adad"]<-mean(a[a$group=="ad-ad",i]=="D")
  prop.amp[i,"adnotad"]<-mean(a[a$group=="ad-notad",i]=="A")
  prop.del[i,"adnotad"]<-mean(a[a$group=="ad-notad",i]=="D")
  }
for(i in sq.genes){
  prop.amp[i,"sqsq"]<-mean(a[a$group=="sq-sq",i]=="A")
  prop.del[i,"sqsq"]<-mean(a[a$group=="sq-sq",i]=="D")
  prop.amp[i,"sqnotsq"]<-mean(a[a$group=="sq-notsq",i]=="A")
  prop.del[i,"sqnotsq"]<-mean(a[a$group=="sq-notsq",i]=="D")
  }

# check NA props align with ad.genes and sq.genes

### plot ###

pdf("cn.pdf",width=3*7)

mymain<-"Proportion of samples with amplification (up) and deletion (down) by LSP group"
x.spec<-c(-0.15,-0.05,0.05,0.15)
myylim<-c(-0.25,0.25)
plot(1:length(genes),xlab="",ylab="Proportion with amp or del",ylim=myylim,xaxt="n",yaxt="n",main=mymain)
for(i in 1:length(genes)){
  mygene<-genes[i]
  if(mygene %in% ad.genes & !mygene %in% sq.genes) {text(i,myylim[1],"ad")}
  if(mygene %in% sq.genes & !mygene %in% ad.genes) {text(i,myylim[1],"sq")}
  segments(x0=i+x.spec[1],y0=0,y1=prop.amp[mygene,"adad"],col="blue",lwd=2)
  segments(x0=i+x.spec[1],y0=0,y1=-prop.del[mygene,"adad"],col="blue",lwd=2)
  segments(x0=i+x.spec[2],y0=0,y1=prop.amp[mygene,"adnotad"],col="green",lwd=2)
  segments(x0=i+x.spec[2],y0=0,y1=-prop.del[mygene,"adnotad"],col="green",lwd=2)

  segments(x0=i+x.spec[3],y0=0,y1=prop.amp[mygene,"sqsq"],col="red",lwd=2)
  segments(x0=i+x.spec[3],y0=0,y1=-prop.del[mygene,"sqsq"],col="red",lwd=2)
  segments(x0=i+x.spec[4],y0=0,y1=prop.amp[mygene,"sqnotsq"],col="orange",lwd=2)
  segments(x0=i+x.spec[4],y0=0,y1=-prop.del[mygene,"sqnotsq"],col="orange",lwd=2)
  }

axis(1,at=1:length(genes),labels=genes,las=2)
axis(2,at=c(myylim[1],0,myylim[2]),labels=abs(c(myylim[1],0,myylim[2])))
abline(h=0,col="gray")
legend("topleft",c("ad-ad","ad-notad","sq-sq","sq-notsq"),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("blue","green","red","orange"),bty="n")

dev.off()

################
### mutation ###
################

# attention to NA--within ad it means no mutation but when I expand the gene list there will be more NA that have different meaning

# 4 groups, plot proportion with Amp, Del, use all genes.  
# Indicate which genes were not in ad or not in sq data.
# -in squam only, or in adeno only (color the gene name).

ad.mut<-data.frame(ad.mut)
sq.mut<-data.frame(sq.mut)
ad.genes<-names(ad.mut)
sq.genes<-names(sq.mut)
genes<-sort(union(names(ad.mut),names(sq.mut)))

ad.mut[,genes[!genes %in% ad.genes]]<-NA
ad.mut<-ad.mut[,genes]
sq.mut[,genes[!genes %in% sq.genes]]<-NA
sq.mut<-sq.mut[,genes]

table(names(ad.mut)==names(sq.mut))
a<-rbind(ad.mut,sq.mut)
a<-data.frame(merge(a,lsp.tcga$l,by=0),row.names=1)

### p-values ###

# ad-notad vs other group (or ad-sq vs other group)?  2x2 table for Amp/not x g1/g2 and table for Del/not x g1/g2
#p.amp<-data.frame(gene=genes,vs.adad=NA,vs.sqsq=NA)
#p.del<-data.frame(gene=genes,vs.adad=NA,vs.sqsq=NA)

### proportions ###

prop.mut<-data.frame(gene=genes,adad=NA,adnotad=NA,sqsq=NA,sqnotsq=NA,row.names=1)
for(i in ad.genes){
  prop.mut[i,"adad"]<-mean(!is.na(a[a$group=="ad-ad",i]))
  prop.mut[i,"adnotad"]<-mean(!is.na(a[a$group=="ad-notad",i]))
  }
for(i in sq.genes){
  prop.mut[i,"sqsq"]<-mean(!is.na(a[a$group=="sq-sq",i]))
  prop.mut[i,"sqnotsq"]<-mean(!is.na(a[a$group=="sq-notsq",i]))
  }

# check NA props align with ad.genes and sq.genes

### plot ###

pdf("mut.pdf",width=3*7)

mymain<-"Proportion of samples with a mutation"
x.spec<-c(-0.15,-0.05,0.05,0.15)
myylim<-c(-.1,0.5)
plot(1:length(genes),xlab="",ylab="Proportion with mutation",ylim=myylim,xaxt="n",yaxt="n",main=mymain)
for(i in 1:length(genes)){
  mygene<-genes[i]
  if(mygene %in% ad.genes & !mygene %in% sq.genes) {text(i,myylim[1],"ad")}
  if(mygene %in% sq.genes & !mygene %in% ad.genes) {text(i,myylim[1],"sq")}
  segments(x0=i+x.spec[1],y0=0,y1=prop.mut[mygene,"adad"],col="blue",lwd=2)
  segments(x0=i+x.spec[2],y0=0,y1=prop.mut[mygene,"adnotad"],col="green",lwd=2)

  segments(x0=i+x.spec[3],y0=0,y1=prop.mut[mygene,"sqsq"],col="red",lwd=2)
  segments(x0=i+x.spec[4],y0=0,y1=prop.mut[mygene,"sqnotsq"],col="orange",lwd=2)
  }

axis(1,at=1:length(genes),labels=genes,las=2)
axis(2,at=c(0,myylim[2]),labels=c(0,myylim[2]))
abline(h=0,col="gray")
legend("topleft",c("ad-ad","ad-notad","sq-sq","sq-notsq"),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("blue","green","red","orange"),bty="n")

dev.off()

### ###