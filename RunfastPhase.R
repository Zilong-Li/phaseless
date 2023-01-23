

########################################
#        sim data
#####################

source("/home/albrecht/tmp2021/2016_02_26/data/anders/misc/ngsAdmix/ngsFatass.R")

set.seed(1)
nSNP<-2000
nInd<-20
dat<-simAdmixture(Nsnp=nSNP,Nind=nInd,c=3,min.maf=0.1,lambda=1/4,region=20)
dat$maf <-apply(dat$geno,2,mean)/2
## simulate sequecing data genotype likelihoods with depth 2
like<-lapply(1:nInd,getLikes2,d=rep(2,nInd),norm=T,geno=dat$geno)



#####################################
## prep data
###############################33
C <- 5
N <- nInd
M <- nSNP
## distance between markers
d<-diff(dat$position)

# random start param
F <- matrix(runif(M*C),ncol=M)
PI <- matrix(runif(M*C),ncol=C)
PI <- PI/rowSums(PI)
par <- list(F=F,PI=PI)
#####################################

 source("../scripts/fastPhase.R")


## run fastPhase
for(i in 1:10){
    cat("iteration",i,"\n")
    logEmission <- emission(like,par$F) #logScale
    trans <- transition(d,par$PI)
    system.time( res <- lapply(1:N,function(x) forwardAndBackwards(x,e=logEmission,trans=trans,M=M,C=C)) )
    cat("log like=", logLike <- sum(sapply(res,function(x) x$totalLikeForwardI)),"\n")
    system.time( par <- updatePar(res,like,C,M,par$F) )
    
}





###############################
# run  on real data
#################################

#FOLDER=/home/albrecht/projects/phaseless/advBinf
#BEAGLE=$FOLDER/prog/beagle.jar
#REF=$FOLDER/exercise2/ref/hg19.fa.gz
#BAMFOLDER=$FOLDER/exercise2/smallBam
#ls $BAMFOLDER/smallNA*.bam > bams.list
#
#angsd -bam bams.list  -SNP_pval 0.001 -doMaf 2 -doMajorMinor 4  -r 1: \
# -doGlf 2  -GL 1 -out results/angsd  -ref $REF
#
#java -jar $BEAGLE like=results/angsd.beagle.gz out=results/imputation
#cp $FOLDER/imputation/hapmap_CEU_r23a_filteredHg19Chr1.tar.gz .
#tar -xf hapmap_CEU_r23a_filteredHg19Chr1.tar.gz



beagle <- read.table("../results/angsd.beagle.gz",head=T,as.is=T)
fed(beagle)
N <- ncol(beagle)/3-1
M <- nrow(beagle)
like <- list()
for(i in 1:N)
    like[[i]] <- t(beagle[,i*3+1:3])

pos <- as.integer( sapply( strsplit(beagle$marker,"_") , function(x) x[2] ) )


C <- 5
d<-diff(pos) / 1e6

# random start param
set.seed(1)
F <- matrix(runif(M*C),ncol=M)
PI <- matrix(runif(M*C),ncol=C)
PI <- PI/rowSums(PI)
par <- list(F=F,PI=PI)
#####################################

 source("../scripts/fastPhase.R")


## run fastPhase
for(i in 1:20){
    cat("iteration",i,"\n")
    logEmission <- emission(like,par$F) #logScale
    trans <- transition(d,par$PI)
    system.time( res <- lapply(1:N,function(x) forwardAndBackwards(x,e=logEmission,trans=trans,M=M,C=C)) )
    cat("log like=", logLike <- sum(sapply(res,function(x) x$totalLikeForwardI)),"\n")
    system.time( par <- updatePar(res,like,C,M,par$F) )  
}

genoCall <- callGeno(res,like,C,M,par$F)

############################3




 
library(snpStats)



################################
#read in the genotype probabilities from beagle
gprobsDat<-read.table("../results/imputation.angsd.beagle.gz.gprobs.gz",head=T,as.is=T)
impuPos<-as.integer(sub("1_","",gprobsDat[,1])) #position

#keep only overlapping positions and convert to a matrix
gprob<-as.matrix(gprobsDat[,-c(1:3)])

#call the genotype with the highest probability
funn<-function(x)
  apply(gprob[,(x-1)*3+1:3],1,which.max)-1

imputa <- t(sapply(1:33, funn))

table(RfastPhase=genoCall$GC,beagle=imputa)

##############################################3
#read in the genotype data. Change the path to where you unpacked the genotypes
pl<-read.plink("../advBinf/hapmap_CEU_r23a_filteredHg19Chr1")
bim<-read.table("../advBinf/hapmap_CEU_r23a_filteredHg19Chr1.bim",as.is=T)
table( keepPL <- bim$V4 %in% impuPos)
indNames<-rownames(pl$genotypes) # individual names in HapMap

indNamesBam<-sapply(strsplit(sub("small","",basename(scan("../bams.list",what="theFck"))),".m"),function(x)x[1]) #individual names in the 33 1000genome
genoHap<-as.integer(pl$genotypes[,keepPL])
genoHap<-matrix(3-genoHap,nrow=60,ncol=sum(keepPL))
rownames(genoHap)<-indNames
genoHap[genoHap==3]<-NA


pos <- bim$V4[keepPL]
mafs<-read.table("../results/angsd.mafs.gz",as.is=T,head=T) #you might have the remove the .gz
mafs<-mafs[mafs$position%in%pos,]
refStrand<-bim[keepPL,5]==mafs$ref
refStrand<-bim[keepPL,5]==mafs$ref
genoHap[rep(refStrand,each=60)]<-2-genoHap[rep(refStrand,each=60)]


#########3
gatk.pos <- read.table("../results/gatk.012.pos")[,2]
sam.pos <-  read.table("../results/sam.012.pos")[,2]
gatk<-read.table("../results/gatk.012")[,-1][gatk.pos %in% pos]
sam<-read.table("../results/sam.012")[,-1][sam.pos %in% pos]





table( keep <- impuPos %in% bim$V4)
dim(genoCall$GC)

#get the mean concordance rate
( concordance <-  c(
bealge=mean(genoHap[indNamesBam,]==imputa[,keep],na.rm=T),
RfastPhase=mean(genoHap[indNamesBam,]==genoCall$GC[,keep],na.rm=T),
samtools=mean(genoHap[indNamesBam,]==sam,na.rm=T),
gatk=mean(genoHap[indNamesBam,]==gatk,na.rm=T)
))
   



