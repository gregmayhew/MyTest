
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

require(xlsx)

#############################
### write 3 squam objects ###
#############################

tmp<-read.xlsx("~/Documents/data/chucks5/tcga/firehose/mutations/from_hawazin/data.file.S7.5.clinical.and.genomic.data.table.xls",sheetIndex=1,startRow=4,header=T,row.names=1)
row.names(tmp)<-gsub("LUSC","TCGA",row.names(tmp))

sq.clin<-tmp[,c(1:5,100)]
names(sq.clin)[names(sq.clin)=="NA."]<-"subtype.tcga"

# copy number in columns 57-99
sq<-tmp[,57:99]
names(sq)<-gsub(".2","",gsub(".1","",names(sq),fixed=T),fixed=T)
#lapply(sq,table,exclude=NULL)
sq.cn<-as.matrix(sq)
sq.cn<-replace(sq.cn,is.na(sq.cn),"N")
sq.cn<-replace(sq.cn,sq.cn=="Amp","A")
sq.cn<-replace(sq.cn,sq.cn=="Del","D")
table(sq.cn,exclude=NULL)

# mutations in columns 13-56
sq<-tmp[,13:56]
sq.mut<-as.matrix(sq)
table(sq.mut,exclude=NULL)

dim(sq.clin); dim(sq.cn); dim(sq.mut)

#############################
### write 3 adeno objects ###
#############################

tmp<-read.csv("~/Documents/data/chucks5/tcga/firehose/mutations/from_hawazin/nature sup 13385-s2 (1).csv",skip=4,header=T,row.names=1,na.strings=c("[Not Available]","[unknown]"),stringsAsFactors=F)
tmp<-tmp[!row.names(tmp)=="",]

ad.clin<-tmp[,c(1:5,36)]

# copy number in columns 68-97
ad<-tmp[,68:97]
names(ad)<-gsub(".2","",gsub(".1","",names(ad),fixed=T),fixed=T)
#lapply(ad,table,exclude=NULL)
ad.cn<-as.matrix(ad)
table(ad.cn,exclude=NULL)

# mutations in columns 12-34
ad<-tmp[,12:34]
ad.mut<-as.matrix(ad)
ad.mut<-replace(ad.mut,ad.mut=="",NA)
ad.mut<-replace(ad.mut,ad.mut=="none",NA)
table(ad.mut,exclude=NULL)

dim(ad.clin); dim(ad.cn); dim(ad.mut)

save(sq.clin,sq.cn,sq.mut,ad.clin,ad.cn,ad.mut,file="paper_downloads.RData")

### ###
