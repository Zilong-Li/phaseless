FATASS2<-function(data,the_ind=c(1),position,maf=NULL,par=NULL,min=0,LD="rsq2",epsilon=0.01,back=0,alim=c(0.01,0.5),optim=1,start=NULL,prune=NULL,ld_adj=TRUE,abs=ifelse(LD=="D",TRUE,FALSE),chr=NULL,diff.state=FALSE,freq,method=NULL,calc.a=F,fix.a=NULL,fix.k2=NULL,moment=NULL){
########################################
if(dim(data)[2]!=length(position))
  stop("number of SNP does not match the number of positions")


if(is.null(maf))
  maf<-apply(data,2,mean,na.rm=TRUE)/2
keep<-!is.na(data[the_ind])&maf>=min&maf<=1-min
maf<-maf[keep]
data<-data[,keep]
ind<-dim(data)[1]
snp<-dim(data)[2]
ind1<-data[the_ind,]
maf<-1-maf
freq<-1-freq
position<-position[keep]
freq<-freq[keep,]
t<-diff(position)
t[t<0]<-1e20
t<-c(0,t,0)

########nyt###############
if(!is.null(chr))
  chr<-chr[keep]
##########################
                                       #udregne LD
mea<-c()
pba<-c()
pBa<-c()
pbA<-c()
pBA<-c()
choose<-c()
##########nyt#######################
if(!is.null(prune)){
  ld<-ld.snp2(data[,snp:1],depth=back)
  keep<-pruning(ld,snp,LD,back,prune,abs)
    maf<-maf[keep]
  data<-data[,keep]
  ind<-dim(data)[1]
  snp<-dim(data)[2]
  ind1<-data[the_ind,]
  if(!is.null(chr))
    chr<-chr[keep]
  position<-position[keep]
  mea<-c()
  cat("pruning done\n")
}

###################################
ep<-error(ind1,epsilon,snp)
S1=emission_Admixture(freq[,1],freq[,2],ep) 

##########
S<-S1

Sk1<-log(S[1,])
Sk2<-log(S[2,])
Sk3<-log(S[3,])
relate_like_Admixture<-function(delta,Sk1,Sk2,Sk3,t,alim=c(0.01,100),fix.a=NULL,fix.k2=NULL,moment=NULL,calc.a=FALSE,phi=0.013){
  
if(calc.a)
  fix.a<-alim[1]
if(!is.null(fix.a))
  delta<-c(fix.a,delta)
if(!is.null(fix.k2))
  delta<-c(delta[1],fix.k2,delta[-1])
if(!is.null(moment))
  delta<-c(delta[1],rev(moment[2:3]))
############################
a<-delta[1]
k0 <- delta[2]
k1 <- delta[3]
k2 <- 1-(k0+k1)
if(calc.a&(k1+2*k2)^2<4*k2){
  cat("qaud error")
  return(1e20)
}
if(sum(c(k0,k2,k1)<0,a<alim[1])>0|sum(c(k0,k1,k2,a-alim[2]+1)>1)>0)
  return(1e20)

if(calc.a)
  a<-calc.a(k0,k1,k2,phi=phi)


IBD01<-(1-exp(-a*t))*k1
IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2

if(k1==1)
  IBD02<-rep(0,length(IBD02))
IBD02[IBD02<1e-15]<-0 #underflow problem loesning
IBD00<-1-IBD01-IBD02
if(k1==0){
IBD10<-IBD01
IBD12<-IBD01
}
else{

IBD10<-IBD01*k0/k1
IBD12<-IBD01*k2/k1
}
IBD11<-1-IBD10-IBD12
IBD21<-IBD01
IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
IBD20[IBD20<1e-15]<-0 #underflow problem loesning

if(k1==1)
  IBD20<-rep(0,length(IBD20))
IBD22<-1-IBD21-IBD20
#apply(cbind(IBD00,IBD10,IBD20,IBD01,IBD11,IBD21,IBD02,IBD12,IBD22),2,min)




  log_IBD0<- log(k0) #start likelihood
  log_IBD1<- log(k1)
  log_IBD2<- log(k2)
  k<-0
  for (i in (2:snp)){
    logp0<-exp(log_IBD0+Sk3[i-1]-k)
    logp1<-exp(log_IBD1+Sk2[i-1]-k)
    logp2<-exp(log_IBD2+Sk1[i-1]-k)
    log_IBD0<-k+log(IBD00[i]*logp0+IBD10[i]*logp1+IBD20[i]*logp2)
    log_IBD1<-k+log(IBD01[i]*logp0+IBD11[i]*logp1+IBD21[i]*logp2)
    log_IBD2<-k+log(IBD02[i]*logp0+IBD12[i]*logp1+IBD22[i]*logp2)
    k<-max(log(logp1)+k,log(logp2)+k,log(logp0)+k)
    if(is.na(k)|is.na(log_IBD0)){
      cat("error in likelihood, delta=\n")
      print(delta)
      break
    }
  }
  logp0 = Sk3[snp] + log_IBD0
  logp1 = Sk2[snp] + log_IBD1
  logp2 = Sk1[snp] + log_IBD2

  k = max(logp0, logp1, logp2)
  l = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k))
#print(c(l,delta))

  return(-l)

}

######
like1<--relate_like_Admixture(c(alim[1],1,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
like2<--relate_like_Admixture(c(alim[1],0,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
min=rep(0,3)
max=rep(1,3)
  conv<-c()
###########nyt##############
#men skal ikke lavet om i c koden
if(is.null(method))
	method<-sample(c("Nelder-Mead","BFGS"),1)
###########################
if(is.null(par)){
  cat("start optim\n")
  conv$value<-Inf
  for(tal in 1:optim){
    temp.start<-diff(c(0,sort(runif(8)),1))[1:8]
    temp.start<-c(1-sum(temp.start[c(1,5)])-sum(temp.start[c(3,6)]),sum(temp.start[c(3,6)]))
    temp.start <- c(runif(1,min=alim[1],max=alim[2]),temp.start)
    if(tal==1&&!is.null(start))
      temp.start=start
###########nyt##############
    if(!is.null(fix.k2))
      temp.start<-temp.start[-2]
    if(!is.null(fix.a)|calc.a)
      temp.start<-temp.start[-1]
    if(!is.null(moment)){
      temp.start<-temp.start[1]
      moment<-moment(maf,1-maf,ind1,ind2)
	}
##################################
    if(length(temp.start)==1)
        conv_temp<-optim(temp.start,relate_like_Admixture,Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim,fix.a=fix.a,fix.k2=fix.k2,moment=moment,method="Brent",calc.a=calc.a,phi=phi,lower=0,upper=1)
    else
        conv_temp<-optim(temp.start,relate_like_Admixture,Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim,fix.a=fix.a,fix.k2=fix.k2,moment=moment,method=method,calc.a=calc.a,phi=phi)
    
  }
  if(conv_temp$value<conv$value)
    conv<-conv_temp
  par<-conv$par
  ##############3nyt
  if(calc.a)
    fix.a<-1  
  if(!is.null(fix.a))
    par<-c(fix.a,par)
  if(!is.null(fix.k2))
    par<-c(par[1],fix.k2,par[-1])
  if(!is.null(moment))
   par<-c(par[1],rev(moment[2:3]))
  if(calc.a)
   par<-c(calc.a(1-sum(par[2:3]),par[3],par[2]),par[-1])  

#####################################
  a <- par[1]#atan(conv$par[1])/pi+0.5
  k<-c(par[2:3],1-sum(par[2:3]))
  k.like<--conv$value
}
else{
a<-par[1]
k<-par[2:4]
k.like<--relate_like_Admixture(c(a,par[2:3]),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
conv<-NULL
}


if(FALSE){
  #conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max)

}
#conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max,Sk=Sk,t=t)
#print(conv)
admix<-k[1]+k[2]/2




obj<-list(k=k,k.like=k.like,admix=admix,a=a,like1=like1,like2=like2,conv=conv,LD=LD,t=t,snp=snp,position=position,S=S,alim=alim,choose=choose,mea=mea,hap=list(pBA,pbA,pBa,pba),S1=S1,maf=maf,chr=chr,conv=conv,freq=freq,geno=ind1)

class(obj)="relateHMM"
obj$decode<-decode(obj)
return(obj)
}



`viterbi_Admixture` <-
function(x){
t<-x$t
Sk<-log(x$S[1:3,])
path<-c()
a <- x$a
k2 <- x$k[3]
k1 <- x$k[2]
k0 <- x$k[1]
snp<-x$snp
#transition
IBD01<-(1-exp(-a*t))*k1
IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2
if(k1==1)
  IBD02<-rep(0,length(IBD02))
IBD02[IBD02<1e-15]<-0 #underflow problem loesning
IBD00<-1-IBD01-IBD02
if(k1==0){
IBD10<-IBD01
IBD12<-IBD01
}
else{
IBD10<-IBD01*k0/k1
IBD12<-IBD01*k2/k1
}
IBD11<-1-IBD10-IBD12
IBD21<-IBD01
IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
IBD20[IBD20<1e-15]<-0 #underflow problem loesning
if(k1==1)
  IBD20<-rep(0,length(IBD20))
IBD22<-1-IBD21-IBD20
#apply(cbind(IBD00,IBD10,IBD20,IBD01,IBD11,IBD21,IBD02,IBD12,IBD22),2,min)

IBD00<-log(IBD00)
IBD10<-log(IBD10)
IBD20<-log(IBD20)
IBD01<-log(IBD01)
IBD11<-log(IBD11)
IBD21<-log(IBD21)
IBD02<-log(IBD02)
IBD12<-log(IBD12)
IBD22<-log(IBD22)

###############################
#viterbi
v0<-max(c(log(k0),log(k1),log(k2)))+Sk[3,1]
v1<-max(c(log(k0),log(k1),log(k2)))+Sk[2,1]
v2<-max(c(log(k0),log(k1),log(k2)))+Sk[1,1]

ptr0<--1
ptr1<--1
ptr2<--1
for (i in (2:snp)){
  v0[i]<-Sk[3,i]+max(c(IBD00[i]+v0[i-1],IBD10[i]+v1[i-1],IBD20[i]+v2[i-1]))
  v1[i]<-Sk[2,i]+max(c(IBD01[i]+v0[i-1],IBD11[i]+v1[i-1],IBD21[i]+v2[i-1]))
  v2[i]<-Sk[1,i]+max(c(IBD02[i]+v0[i-1],IBD12[i]+v1[i-1],IBD22[i]+v2[i-1]))
  ptr0[i]<-which.max(c(v0[i-1]+IBD00[i],v1[i-1]+IBD10[i],v2[i-1]+IBD20[i]))-1
  ptr1[i]<-which.max(c(v0[i-1]+IBD01[i],v1[i-1]+IBD11[i],v2[i-1]+IBD21[i]))-1
  ptr2[i]<-which.max(c(v0[i-1]+IBD02[i],v1[i-1]+IBD12[i],v2[i-1]+IBD22[i]))-1
}

  l = max(c(v0[snp],v1[snp],v2[snp]))
#print(l      )
pi<-rep(0,snp)
pi[snp]<-which.max(c(v0[snp],v1[snp],v2[snp]))-1

for(i in snp:2){
  if(pi[i]==2)
    pi[i-1]<-ptr2[i]
  else if(pi[i]==1)
    pi[i-1]<-ptr1[i]
  else
    pi[i-1]<-ptr0[i]
} 
return(pi)
}






print.relateHMM<-function(x){
 if (!inherits(x, "relateHMM")) 
        stop("Not an object of class relate!")
print(structure(list(k=round(x$k,2),k.like=round(x$k.like,3),admix=signif(x$admix,3),a=signif(x$a,3),like1=round(x$like1,1),like2=round(x$like2,1)),class="power.htest"))
}


plot.relateHMM<-function(x,col=1:3,lwd=2,chr=NULL,...){
post<-x$decode$post
path<-x$decode$path
pos<-x$position
if(is.null(x$chr)){
  plot(pos[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
  points(pos[1:(x$snp-1)],rep(1,x$snp-1),col=3-path[-(x$snp-1)],pch="|")
}
else{
  if(!is.null(chr)){
    pos<-pos[x$chr==chr]
    post<-post[x$chr==chr,]
    plot(pos[-length(pos)],post[-1,1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
    lines(pos[-length(pos)],post[-1,2],col=col[2],lwd=lwd)
    lines(pos[-length(pos)],post[-1,3],col=col[3],lwd=lwd)
    
  }
  else{
    C<-names(table(x$chr))
    m<-c(0,cumsum(tapply(pos,x$chr,max)))
    pos2<-rep(NA,length(pos))
    for(tal in 1:length(C))
      pos2[x$chr==C[tal]]<-pos[x$chr==C[tal]]+m[tal]
    plot(pos2[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
    lines(pos2[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
    lines(pos2[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
    abline(v=m)
    for(tal in 1:length(C))
      text(m[tal]+diff(m)[tal]/2,0.5,C[tal],col="gray")
  }
}
invisible(cbind(post,pos))
}


decode<-function(x){
t<-x$t
Sk<-log(x$S[1:3,])
post<-c()
path<-c()
a <- x$a
k0 <- x$k[1]
k1 <- x$k[2]
k2 <- x$k[3]
snp<-x$snp
#transition
IBD01<-(1-exp(-a*t))*k1
IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2
if(k1==1)
  IBD02<-rep(0,length(IBD02))
IBD02[IBD02<1e-15]<-0 #underflow problem loesning
IBD00<-1-IBD01-IBD02
if(k1==0){
IBD10<-IBD01
IBD12<-IBD01
}
else{
IBD10<-IBD01*k0/k1
IBD12<-IBD01*k2/k1
}
IBD11<-1-IBD10-IBD12
IBD21<-IBD01
IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
IBD20[IBD20<1e-15]<-0 #underflow problem loesning
if(k1==1)
  IBD20<-rep(0,length(IBD20))
IBD22<-1-IBD21-IBD20
#apply(cbind(IBD00,IBD10,IBD20,IBD01,IBD11,IBD21,IBD02,IBD12,IBD22),2,min)


###############################
 #start forward
  log_IBD0<-c()
  log_IBD1<-c()
  log_IBD2<-c()
  log_IBD0[1]<- log(k0)+Sk[3,1]
  log_IBD1[1]<- log(k1)+Sk[2,1]
  log_IBD2[1]<- log(k2)+Sk[1,1]
  k<-0
  for (i in (2:snp)){
    logp0<-exp(log_IBD0[i-1]-k)
    logp1<-exp(log_IBD1[i-1]-k)
    logp2<-exp(log_IBD2[i-1]-k)
    log_IBD0[i]<-k+log(IBD00[i]*logp0+IBD10[i]*logp1+IBD20[i]*logp2)+Sk[3,i]
    log_IBD1[i]<-k+log(IBD01[i]*logp0+IBD11[i]*logp1+IBD21[i]*logp2)+Sk[2,i]
    log_IBD2[i]<-k+log(IBD02[i]*logp0+IBD12[i]*logp1+IBD22[i]*logp2)+Sk[1,i]
    k<-max(log_IBD0[i],log_IBD1[i],log_IBD2[i])
   }
  logp0 = log_IBD0[snp]
  logp1 = log_IBD1[snp]
  logp2 = log_IBD2[snp]
  k = max(logp0, logp1, logp2)
  l = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k))
##########################

  bk_log_IBD0<-c()
  bk_log_IBD1<-c()
  bk_log_IBD2<-c()
  bk_log_IBD0[snp]<-0
  bk_log_IBD1[snp]<-0
  bk_log_IBD2[snp]<-0
  k<-0
  for (i in ((snp-1):1)){
    p0<-exp(Sk[3,i+1]+bk_log_IBD0[i+1]-k)
    p1<-exp(Sk[2,i+1]+bk_log_IBD1[i+1]-k)
    p2<-exp(Sk[1,i+1]+bk_log_IBD2[i+1]-k)
    bk_log_IBD0[i]<-k+log(p2*IBD02[i+1]+p1*IBD01[i+1]+p0*IBD00[i+1])
    bk_log_IBD1[i]<-k+log(p2*IBD12[i+1]+p1*IBD11[i+1]+p0*IBD10[i+1])
    bk_log_IBD2[i]<-k+log(p2*IBD22[i+1]+p1*IBD21[i+1]+p0*IBD20[i+1])
    k<-max(bk_log_IBD2[i],bk_log_IBD1[i],bk_log_IBD0[i])
   }
  logp0 =  bk_log_IBD0[1]+ Sk[3,1]+log(k0)
  logp1 =  bk_log_IBD1[1]+ Sk[2,1]+log(k1)
  logp2 =  bk_log_IBD2[1]+ Sk[1,1]+log(k2)
  k = max(logp0, logp1, logp2)
 # l2 = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k))
# c(l,l2)
###########################
post<-cbind(exp(log_IBD2+bk_log_IBD2-l),exp(log_IBD1+bk_log_IBD1-l),exp(log_IBD0+bk_log_IBD0-l))
apply(post,1,sum)

path<-viterbi_Admixture(x)
#####
obj=list(path=path,post=post)
class(obj)="decode"

return(obj)
}









error<-function(geno,epsilon,snp){
  ep<-matrix(NA,ncol=3,nrow=snp)

ep[geno==0,1]<-(1-epsilon)^2
ep[geno==1,1]<-(1-epsilon)*epsilon
ep[geno==2,1]<-(epsilon)^2

ep[geno==0,2]<-2*(1-epsilon)*epsilon
ep[geno==1,2]<-(1-epsilon)^2+epsilon^2
ep[geno==2,2]<-2*(1-epsilon)*epsilon

ep[geno==0,3]<-(epsilon)^2
ep[geno==1,3]<-(1-epsilon)*epsilon
ep[geno==2,3]<-(1-epsilon)^2

return(ep)
}

simAdmixtureFreq<-function(Nsnp=1000,shape=c(2,2),c=10){
prob=rbeta(Nsnp,shape1=shape[1],shape2=shape[2])
prob2=rbeta(Nsnp,shape1=prob*(c-1),shape2=c+prob-c*prob-1)
return(cbind(prob,prob2))
}


AdmixtureState<-function(Nsnp=1000,Nind=20,admixture=runif(Nind*2),lambda=50){
state<-matrix(NA,ncol=Nsnp,nrow=Nind*2)
for(tal in 1:(Nind*2)){
p<-sample(1:2,Nsnp/10,prob=c(admixture[tal],1-admixture[tal]),replace=T)
l<-0
while(sum(l)<Nsnp)
l<-rgeom(Nsnp/10,lambda)*100
state[tal,]<-rep(p,l)[1:Nsnp]
}
if(any(is.na(state)))
   stop("fucked shit")
return(state)
}

AdmixtureGeno<-function(state,freq){
geno1<-matrix(rbinom(prod(dim(state)),1,prob=rep(freq[,1],each=nrow(state))),ncol=ncol(state))
geno2<-matrix(rbinom(prod(dim(state)),1,prob=rep(freq[,2],each=nrow(state))),ncol=ncol(state))
geno1[state==2]<-geno2[state==2]
geno1<-geno1[1:(nrow(state)/2)*2-1,]+geno1[1:(nrow(state)/2)*2,]
return(geno1)
}

emission_Admixture<-function(freq1,freq2,ep){
  S<-cbind(freq1^2,freq1*freq2,freq2^2)*ep[,1] #AA
  S<-S+cbind(2*freq1*(1-freq1),freq1*(1-freq2)+freq2*(1-freq1),2*freq2*(1-freq2))*ep[,2]#Aa
  S<-S+cbind((1-freq1)^2,(1-freq1)*(1-freq2),(1-freq2)^2)*ep[,3]#aa
  return(t(S))
}


#data<-geno
#min<-0.01
#the_ind=1
#par=NULL;min=0;LD="rsq2";epsilon=0.01;back=0;alim=c(0.01,1);optim=1;start=NULL;prune=NULL;ld_adj=TRUE;abs=ifelse(LD=="D",TRUE,FALSE);chr=NULL;diff.state=FALSE;calc.a=F;fix.a=NULL;fix.k2=NULL;moment=NULL



`FATASS` <-
function(data,the_ind=c(1),position,maf=NULL,par=NULL,min=0,LD="rsq2",epsilon=0.01,back=0,alim=c(0.01,0.5),optim=1,start=NULL,prune=NULL,ld_adj=TRUE,abs=ifelse(LD=="D",TRUE,FALSE),chr=NULL,diff.state=FALSE,freq,method=NULL,calc.a=F,fix.a=NULL,fix.k2=NULL,moment=NULL){
########################################
if(dim(data)[2]!=length(position))
  stop("number of SNP does not match the number of positions")


if(is.null(maf))
  maf<-apply(data,2,mean,na.rm=TRUE)/2
keep<-!is.na(data[the_ind])&maf>=min&maf<=1-min
maf<-maf[keep]
data<-data[,keep]
ind<-dim(data)[1]
snp<-dim(data)[2]
ind1<-data[the_ind,]
maf<-1-maf
freq<-1-freq
position<-position[keep]
freq<-freq[keep,]
t<-diff(position)
t[t<0]<-1e20
t<-c(0,t,0)

########nyt###############
if(!is.null(chr))
  chr<-chr[keep]
##########################
                                       #udregne LD
mea<-c()
pba<-c()
pBa<-c()
pbA<-c()
pBA<-c()
choose<-c()
##########nyt#######################
if(!is.null(prune)){
  ld<-ld.snp2(data[,snp:1],depth=back)
  keep<-pruning(ld,snp,LD,back,prune,abs)
    maf<-maf[keep]
  data<-data[,keep]
  ind<-dim(data)[1]
  snp<-dim(data)[2]
  ind1<-data[the_ind,]
  if(!is.null(chr))
    chr<-chr[keep]
  position<-position[keep]
  mea<-c()
  cat("pruning done\n")
}

###################################
ep<-error(ind1,epsilon,snp)
S1=emission_Admixture(freq[,1],freq[,2],ep) 

##########
S<-S1

Sk1<-log(S[1,])
Sk2<-log(S[2,])
Sk3<-log(S[3,])
relate_like_Admixture<-function(delta,Sk1,Sk2,Sk3,t,alim=c(0.01,100),fix.a=NULL,fix.k2=NULL,moment=NULL,calc.a=FALSE,phi=0.013){
  
if(calc.a)
  fix.a<-alim[1]
if(!is.null(fix.a))
  delta<-c(fix.a,delta)
if(!is.null(fix.k2))
  delta<-c(delta[1],fix.k2,delta[-1])
if(!is.null(moment))
  delta<-c(delta[1],rev(moment[2:3]))
############################
a<-delta[1]
k0 <- delta[2]
k1 <- delta[3]
k2 <- 1-(k0+k1)
if(calc.a&(k1+2*k2)^2<4*k2){
  cat("qaud error")
  return(1e20)
}
if(sum(c(k0,k2,k1)<0,a<alim[1])>0|sum(c(k0,k1,k2,a-alim[2]+1)>1)>0)
  return(1e20)

if(calc.a)
  a<-calc.a(k0,k1,k2,phi=phi)


IBD01<-(1-exp(-a*t))*k1
IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2

if(k1==1)
  IBD02<-rep(0,length(IBD02))
IBD02[IBD02<1e-15]<-0 #underflow problem loesning
IBD00<-1-IBD01-IBD02
if(k1==0){
IBD10<-IBD01
IBD12<-IBD01
}
else{

IBD10<-IBD01*k0/k1
IBD12<-IBD01*k2/k1
}
IBD11<-1-IBD10-IBD12
IBD21<-IBD01
IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
IBD20[IBD20<1e-15]<-0 #underflow problem loesning

if(k1==1)
  IBD20<-rep(0,length(IBD20))
IBD22<-1-IBD21-IBD20
#apply(cbind(IBD00,IBD10,IBD20,IBD01,IBD11,IBD21,IBD02,IBD12,IBD22),2,min)




  log_IBD0<- log(k0) #start likelihood
  log_IBD1<- log(k1)
  log_IBD2<- log(k2)
  k<-0
  for (i in (2:snp)){
    logp0<-exp(log_IBD0+Sk3[i-1]-k)
    logp1<-exp(log_IBD1+Sk2[i-1]-k)
    logp2<-exp(log_IBD2+Sk1[i-1]-k)
    log_IBD0<-k+log(IBD00[i]*logp0+IBD10[i]*logp1+IBD20[i]*logp2)
    log_IBD1<-k+log(IBD01[i]*logp0+IBD11[i]*logp1+IBD21[i]*logp2)
    log_IBD2<-k+log(IBD02[i]*logp0+IBD12[i]*logp1+IBD22[i]*logp2)
    k<-max(log(logp1)+k,log(logp2)+k,log(logp0)+k)
    if(is.na(k)|is.na(log_IBD0)){
      cat("error in likelihood, delta=\n")
      print(delta)
      break
    }
  }
  logp0 = Sk3[snp] + log_IBD0
  logp1 = Sk2[snp] + log_IBD1
  logp2 = Sk1[snp] + log_IBD2

  k = max(logp0, logp1, logp2)
  l = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k))
#print(c(l,delta))

  return(-l)

}

######
like1<--relate_like_Admixture(c(alim[1],1,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
like2<--relate_like_Admixture(c(alim[1],0,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
min=rep(0,3)
max=rep(1,3)
  conv<-c()
###########nyt##############
#men skal ikke lavet om i c koden
if(is.null(method))
	method<-sample(c("Nelder-Mead","BFGS"),1)
###########################
if(is.null(par)){
  cat("start optim\n")
  conv$value<-Inf
  for(tal in 1:optim){
    temp.start<-diff(c(0,sort(runif(8)),1))[1:8]
    temp.start<-c(1-sum(temp.start[c(1,5)])-sum(temp.start[c(3,6)]),sum(temp.start[c(3,6)]))
    temp.start <- c(runif(1,min=alim[1],max=alim[2]),temp.start)
    if(tal==1&&!is.null(start))
      temp.start=start
###########nyt##############
    if(!is.null(fix.k2))
      temp.start<-temp.start[-2]
    if(!is.null(fix.a)|calc.a)
      temp.start<-temp.start[-1]
    if(!is.null(moment)){
      temp.start<-temp.start[1]
moment<-moment(maf,1-maf,ind1,ind2)
	}
##################################
    conv_temp<-optim(temp.start,relate_like_Admixture,Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim,fix.a=fix.a,fix.k2=fix.k2,moment=moment,method=method,calc.a=calc.a,phi=phi)
  }
  if(conv_temp$value<conv$value)
    conv<-conv_temp
  par<-conv$par
  ##############3nyt
  if(calc.a)
    fix.a<-1  
  if(!is.null(fix.a))
    par<-c(fix.a,par)
  if(!is.null(fix.k2))
    par<-c(par[1],fix.k2,par[-1])
  if(!is.null(moment))
   par<-c(par[1],rev(moment[2:3]))
  if(calc.a)
   par<-c(calc.a(1-sum(par[2:3]),par[3],par[2]),par[-1])  

#####################################
  a <- par[1]#atan(conv$par[1])/pi+0.5
  k<-c(par[2:3],1-sum(par[2:3]))
  k.like<--conv$value
}
else{
a<-par[1]
k<-par[2:4]
k.like<--relate_like_Admixture(c(a,par[2:3]),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
conv<-NULL
}


if(FALSE){
  #conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max)

}
#conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max,Sk=Sk,t=t)
#print(conv)
admix<-k[1]+k[2]/2




obj<-list(k=k,k.like=k.like,admix=admix,a=a,like1=like1,like2=like2,conv=conv,LD=LD,t=t,snp=snp,position=position,S=S,alim=alim,choose=choose,mea=mea,hap=list(pBA,pbA,pBa,pba),S1=S1,maf=maf,chr=chr,conv=conv,freq=freq,geno=ind1)

class(obj)="relateHMM"
obj$decode<-decode(obj)
return(obj)
}

simAdmixture<-function(Nsnp=5000,Nind=20,c=5,min.maf=0.1,lambda=1/9,region=10){

ft<-simAdmixtureFreq(Nsnp*2,shape=c(2,2),c=c)
sam<-sample(c(TRUE,FALSE),Nsnp*2,replace=T)
f<-ft
f[sam,]<-ft[sam,]
f[!sam,]<-ft[!sam,2:1]
f<-f[apply(abs(f-0.5),1,min)< 0.5-min.maf,]
f<-f[1:Nsnp,]
#plot(f)

st<-AdmixtureState(lam=lambda,Nsnp=Nsnp)
geno<-AdmixtureGeno(st,f)
state<-st[1:Nind*2-1,]+st[1:Nind*2,]-1
# plot(position,rep(0,length(position)),col=state[1,],pch="|")
position<-1:Nsnp/Nsnp*region

return(list(geno=geno,position=position,state=state,st=st,f=f))
}


########################################
if(FALSE){

source("admixure.R")
set.seed(1)
dat<-simAdmixture(Nsnp=5000,Nind=20,c=5,min.maf=0.1,lambda=1/9)

A<-FATASS(dat$geno,1,pos=dat$position,freq=dat$f,alim=c(0.0005,0.50),min=0,epsilon=0.01);A;plot(A)
points(dat$position,rep(0,length(dat$position)),col=dat$state[1,],pch="|")



#cbind(t(A$S[,1:10]),geno[1,1:10],f[1:10,],state[1,1:10])







}





`FATASSngs` <-
function(llike,the_ind=c(1),position,maf=NULL,par=NULL,min=0,LD="rsq2",epsilon=0.01,back=0,alim=c(0.01,0.5),optim=1,start=NULL,prune=NULL,ld_adj=TRUE,abs=ifelse(LD=="D",TRUE,FALSE),chr=NULL,diff.state=FALSE,freq,method=NULL,calc.a=F,fix.a=NULL,fix.k2=NULL,moment=NULL){
#the_ind=c(1);maf=NULL;par=NULL;min=0;LD="rsq2";epsilon=0.01;back=0;alim=c(0.01,0.5);optim=1;start=NULL;prune=NULL;ld_adj=TRUE;abs=ifelse(LD=="D",TRUE,FALSE);chr=NULL;diff.state=FALSE;method=NULL;calc.a=F;fix.a=NULL;fix.k2=NULL;moment=NULL
########################################

snp<-nrow(llike)/3
ind<-ncol(llike)
if(is.null(maf))
  maf<-apply(data,2,mean,na.rm=TRUE)/2
keep<-rep(TRUE,)
maf<-maf[keep]
llike<-llike[rep(keep,each=3),]
snp<-nrow(llike)/3
ind<-ncol(llike)
ind1<-llike[,the_ind]
maf<-1-maf
freq<-1-freq
position<-position[keep]
freq<-freq[keep,]
t<-diff(position)
t[t<0]<-1e20
t<-c(0,t,0)

########nyt###############
if(!is.null(chr))
  chr<-chr[keep]
##########################
                                       #udregne LD
mea<-c()
pba<-c()
pBa<-c()
pbA<-c()
pBA<-c()
choose<-c()
##########nyt#######################
if(!is.null(prune)){
  ld<-ld.snp2(data[,snp:1],depth=back)
  keep<-pruning(ld,snp,LD,back,prune,abs)
    maf<-maf[keep]
  data<-data[,keep]
  ind<-dim(data)[1]
  snp<-dim(data)[2]
  ind1<-data[the_ind,]
  if(!is.null(chr))
    chr<-chr[keep]
  position<-position[keep]
  mea<-c()
  cat("pruning done\n")
}

###################################
 ep<-matrix(llike[,the_ind],ncol=3,by=T)
S1=emission_Admixture(freq[,1],freq[,2],ep) 

##########
S<-S1

Sk1<-log(S[1,])
Sk2<-log(S[2,])
Sk3<-log(S[3,])
relate_like_Admixture<-function(delta,Sk1,Sk2,Sk3,t,alim=c(0.01,100),fix.a=NULL,fix.k2=NULL,moment=NULL,calc.a=FALSE,phi=0.013){
  
if(calc.a)
  fix.a<-alim[1]
if(!is.null(fix.a))
  delta<-c(fix.a,delta)
if(!is.null(fix.k2))
  delta<-c(delta[1],fix.k2,delta[-1])
if(!is.null(moment))
  delta<-c(delta[1],rev(moment[2:3]))
############################
a<-delta[1]
k0 <- delta[2]
k1 <- delta[3]
k2 <- 1-(k0+k1)
if(calc.a&(k1+2*k2)^2<4*k2){
  cat("qaud error")
  return(1e20)
}
if(sum(c(k0,k2,k1)<0,a<alim[1])>0|sum(c(k0,k1,k2,a-alim[2]+1)>1)>0)
  return(1e20)

if(calc.a)
  a<-calc.a(k0,k1,k2,phi=phi)


IBD01<-(1-exp(-a*t))*k1
IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2

if(k1==1)
  IBD02<-rep(0,length(IBD02))
IBD02[IBD02<1e-15]<-0 #underflow problem loesning
IBD00<-1-IBD01-IBD02
if(k1==0){
IBD10<-IBD01
IBD12<-IBD01
}
else{

IBD10<-IBD01*k0/k1
IBD12<-IBD01*k2/k1
}
IBD11<-1-IBD10-IBD12
IBD21<-IBD01
IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
IBD20[IBD20<1e-15]<-0 #underflow problem loesning

if(k1==1)
  IBD20<-rep(0,length(IBD20))
IBD22<-1-IBD21-IBD20
#apply(cbind(IBD00,IBD10,IBD20,IBD01,IBD11,IBD21,IBD02,IBD12,IBD22),2,min)




  log_IBD0<- log(k0) #start likelihood
  log_IBD1<- log(k1)
  log_IBD2<- log(k2)
  k<-0
  for (i in (2:snp)){
    logp0<-exp(log_IBD0+Sk3[i-1]-k)
    logp1<-exp(log_IBD1+Sk2[i-1]-k)
    logp2<-exp(log_IBD2+Sk1[i-1]-k)
    log_IBD0<-k+log(IBD00[i]*logp0+IBD10[i]*logp1+IBD20[i]*logp2)
    log_IBD1<-k+log(IBD01[i]*logp0+IBD11[i]*logp1+IBD21[i]*logp2)
    log_IBD2<-k+log(IBD02[i]*logp0+IBD12[i]*logp1+IBD22[i]*logp2)
    k<-max(log(logp1)+k,log(logp2)+k,log(logp0)+k)
    if(is.na(k)|is.na(log_IBD0)){
      cat("error in likelihood, delta=\n")
      print(delta)
      break
    }
  }
  logp0 = Sk3[snp] + log_IBD0
  logp1 = Sk2[snp] + log_IBD1
  logp2 = Sk1[snp] + log_IBD2

  k = max(logp0, logp1, logp2)
  l = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k))
#print(c(l,delta))

  return(-l)

}

######
like1<--relate_like_Admixture(c(alim[1],1,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
like2<--relate_like_Admixture(c(alim[1],0,0),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
min=rep(0,3)
max=rep(1,3)
  conv<-c()
###########nyt##############
#men skal ikke lavet om i c koden
if(is.null(method))
	method<-sample(c("Nelder-Mead","BFGS"),1)
###########################
if(is.null(par)){
  cat("start optim\n")
  conv$value<-Inf
  for(tal in 1:optim){
    temp.start<-diff(c(0,sort(runif(8)),1))[1:8]
    temp.start<-c(1-sum(temp.start[c(1,5)])-sum(temp.start[c(3,6)]),sum(temp.start[c(3,6)]))
    temp.start <- c(runif(1,min=alim[1],max=alim[2]),temp.start)
    if(tal==1&&!is.null(start))
      temp.start=start
###########nyt##############
    if(!is.null(fix.k2))
      temp.start<-temp.start[-2]
    if(!is.null(fix.a)|calc.a)
      temp.start<-temp.start[-1]
    if(!is.null(moment)){
      temp.start<-temp.start[1]
moment<-moment(maf,1-maf,ind1,ind2)
	}
##################################
    conv_temp<-optim(temp.start,relate_like_Admixture,Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim,fix.a=fix.a,fix.k2=fix.k2,moment=moment,method=method,calc.a=calc.a,phi=phi)
  }
  if(conv_temp$value<conv$value)
    conv<-conv_temp
  par<-conv$par
  ##############3nyt
  if(calc.a)
    fix.a<-1  
  if(!is.null(fix.a))
    par<-c(fix.a,par)
  if(!is.null(fix.k2))
    par<-c(par[1],fix.k2,par[-1])
  if(!is.null(moment))
   par<-c(par[1],rev(moment[2:3]))
  if(calc.a)
   par<-c(calc.a(1-sum(par[2:3]),par[3],par[2]),par[-1])  

#####################################
  a <- par[1]#atan(conv$par[1])/pi+0.5
  k<-c(par[2:3],1-sum(par[2:3]))
  k.like<--conv$value
}
else{
a<-par[1]
k<-par[2:4]
k.like<--relate_like_Admixture(c(a,par[2:3]),Sk1=Sk1,Sk2=Sk2,Sk3=Sk3,t=t,alim=alim)
conv<-NULL
}


if(FALSE){
  #conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max)

}
#conv<-nlminb(temp.start,relate_like_HMM_k,lower=min,upper=max,Sk=Sk,t=t)
#print(conv)
admix<-k[1]+k[2]/2




obj<-list(k=k,k.like=k.like,admix=admix,a=a,like1=like1,like2=like2,conv=conv,LD=LD,t=t,snp=snp,position=position,S=S,alim=alim,choose=choose,mea=mea,hap=list(pBA,pbA,pBa,pba),S1=S1,maf=maf,chr=chr,conv=conv,freq=freq,geno=ind1)

class(obj)="relateHMM"
obj$decode<-decode(obj)
return(obj)
}




################################################
#   stuff for simulating sequencing data and genotype likelihoods
################################################
plotter<-function(x,col=1:8,lwd=2,...){
  post<-x$decode$post
  path<-x$decode$path
  pos<-x$position

  plot(pos[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
  legend(6,0.5,paste("ancestral state",0:2),col=1:3,lty=1)

}



genoCall<-function(x,cut=0){
   ind<-ncol(x)

   for(i in 1:ind){
     w<-apply(matrix(x[,i],3),2,which.max)
     mat<-matrix(0,nrow=3,ncol=length(w))
     mat[(1:length(w)-1)*3+w]<-1
     pp<-x[(1:length(w)-1)*3+w,i]
     mat[,pp<cut]<-1
     x[,i]<-mat
   }
   x
}

################################3


genoCall<-function(x,cut=0){
   ind<-ncol(x)

   for(i in 1:ind){
     w<-apply(matrix(x[,i],3),2,which.max)
     mat<-matrix(0,nrow=3,ncol=length(w))
     mat[(1:length(w)-1)*3+w]<-1
     pp<-x[(1:length(w)-1)*3+w,i]
     mat[,pp<cut]<-1
     x[,i]<-mat
   }
   x
}
writeBeagleFormat<-function(x,filename){
  options(scipen=999)
  nInd<-length(x)
  if(class(x[[1]])=="matrix"){
    nSNP<-ncol(x[[1]])
    mat<-matrix(NA,nrow=nSNP,ncol=3*nInd)
    for(l in 1:nInd-1)
      mat[,l*3+1:3]<-t(x[[l+1]])
  }
  else{#if genotypes
    nSNP<-length(x[[1]])
    mat<-matrix(NA,nrow=nSNP,ncol=3*nInd)
    for(l in 1:nInd-1){
      minimat<-matrix(0,nrow=3,ncol=nSNP)
      minimat[(2-x[[l+1]])+(1:nSNP-1)*3+1]<-1
      mat[,l*3+1:3]<-t(minimat)

    }
  }
  header<-c("marker","allele1","allele2",paste("ind",rep(1:nInd,each=3),sep=""))
  res<-cbind(paste("SNP",1:nSNP,sep=""),"0","1",mat) 
  write.table(rbind(header),file=filename,row=FALSE,col=FALSE,qu=FALSE,sep="\t")
  write.table(res,file=filename,row=FALSE,col=FALSE,qu=FALSE,sep="\t",append=T)
}

simFreq<-function(x,nSNP,fac=2,minMaf=0.05,maxMaf=0.5){
  npop<-nrow(x)
  nInd<-ncol(x)
  f<-runif(nSNP,min=minMaf,max=maxMaf)
  ffun<-function(x){
    f<-rnorm(nSNP,mean=f,sd=(f*(1-f))^2*fac)
    f<-ifelse(f<0|f>1,0.3,f)
  }
 sapply(1:npop,ffun)
}

getLikes<-function(x,d=5,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}

getLikes2<-function(y,dep,geno,e=0.01,norm=FALSE){
  x<-geno[y,]
  d<-dep[y]
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}


simGeno<-function(x,f){
  len<-nrow(f)
  F<-f%*%x
  rbinom(len,2,F)
}

loglikeKnown<-function(q,f,geno,npop=2){
  F<-1-t(f)
  loglikeInd<-function(g,qPart){
    #sum(g*log(qPart%*%f)+(2-g)*log(qPart%*%(1-f)))
    fPart<-F%*%qPart
    sum((2-g)*log(fPart)+(g)*log(1-fPart))+log(2)*sum(g==1)
  }
  loglike<-0
  for(tal in 1:length(geno))
    loglike<-loglike+loglikeInd(geno[[tal]],q[,tal])
 # loglike<-sum(unlist(mclapply(1:length(geno),function(tal) loglikeInd(geno[[tal]],q[,tal]),mc.cores=10)))
  return(loglike)
}






loglikeUnKnown<-function(q,f,like,npop=2){
  F<-1-t(f)

  loglikeIndlike<-function(lPart,qPart){
    fPart<-F%*%qPart
    sum(log((fPart^2*lPart[1,]+2*fPart*(1-fPart)*lPart[2,]+(1-fPart)^2*lPart[3,])))
  }
  loglike<-0
  for(tal in 1:length(like))
    loglike<-loglike+loglikeIndlike(like[[tal]],q[,tal])
 # loglike<-sum(unlist(mclapply(1:length(geno),function(tal) loglikeInd(geno[[tal]],q[,tal]),mc.cores=10)))
  return(loglike)
}



loglikeUnKnown3<-function(q,f,llike,npop=2){
  F<-1-t(f)

  loglikeIndlike<-function(lPart,qPart){
    fPart<-F%*%qPart
    sum(log((fPart^2*lPart[1,]+2*fPart*(1-fPart)*lPart[2,]+(1-fPart)^2*lPart[3,])))
  }
  loglike<-0
  for(tal in 1:ncol(llike))
    loglike<-loglike+loglikeIndlike(matrix(llike[,tal],3),q[,tal])
 # loglike<-sum(unlist(mclapply(1:length(geno),function(tal) loglikeInd(geno[[tal]],q[,tal]),mc.cores=10)))
  return(loglike)
}


kk<-c()
loglikeUnKnownWrap<-function(x,like,npop=2){
  kk<<-x
  ind<-length(like)
  sites<-length(like[[1]])
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
    f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  if(any(x<0)|any(x>1)|any(q>1))
    return(Inf)


  ll<--loglikeUnKnown(q=q,f=f,like=like,npop=npop)
  cat(ll,"\n")
  ll
}

loglikeUnKnownWrap3<-function(x,llike,npop=2){
  ind<-ncol(llike)
  sites<-nrow(llike)/3
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  if(any(x<0)|any(x>1)|any(q>1)|any(q<0))
    return(1e20)


  ll<--loglikeUnKnown3(q=q,f=f,llike=llike,npop=npop)
 # cat(ll,"\n")
  ll
}

loglikeUnKnownWrap4<-function(x,llike,npop=2){
  ind<-ncol(llike)
  sites<-nrow(llike)/3
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  if(any(x<0)|any(x>1)|any(q>1)|any(q<0))
    return(1e20)


logLike<-0
  
   res<-likeCpp(f=1-f,l=llike,K=npop,q=q,ind,sites,logLike)
-res$logLike
}






kk<-c()
loglikeKnownWrap<-function(x,geno,npop=2){
  kk<<-x
  ind<-length(geno)
  sites<-length(geno[[1]])
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  if(any(x<0)|any(x>1)|any(q>1)){
    cat("inf\n")
    return(Inf)
  }


  ll<--loglikeKnown(q=q,f=f,geno=geno,npop=npop)
#   cat(ll,"\n")
  ll
}


writeGeno<-function(geno,fileName="file.geno"){

 genoMat<-matrix(unlist(geno),ncol=length(geno))
 write.table(genoMat,file=fileName,sep="",col=FALSE,qu=FALSE,row=FALSE)

}
em<-function(x,geno,npop=2){
  ind<-length(geno)
  sites<-length(geno[[1]])
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  genoMat<-matrix(unlist(geno),ncol=length(geno))
 

  aNorm<-matrix(NA,nrow=ind,ncol=sites)
  bNorm<-matrix(NA,nrow=ind,ncol=sites)
  for(tal in 1:sites){
    aNorm[,tal]<-t(q)%*%f[,tal]
    bNorm[,tal]<-t(q)%*%(1-f[,tal])
  }

  updateF<-function(j,k){
    ag<-genoMat[j,]*q[k,]*f[k,j]/aNorm[,j]
    bg<-(2-genoMat[j,])*q[k,]*(1-f[k,j])/bNorm[,j]
    fnew<-sum(ag)/(sum(ag)+sum(bg))
    fnew
  }
##  (f[1,1]<-updateF(1,1))
  updateQ<-function(i,k){
    ag<-genoMat[,i]*q[k,i]*f[k,]/aNorm[i,]
    bg<-(2-genoMat[,i])*q[k,i]*(1-f[k,])/bNorm[i,]
    qnew<-sum(bg+ag)/(2*sites)
    qnew
  }

  
  fnew<-f
  qnew<-q
  for(k in 1:npop)
    for(j in 1:sites)
      fnew[k,j]<-updateF(j,k)

  
  for(k in 1:(npop-1))
    for(i in 1:ind)
      qnew[k,i]<-updateQ(i,k)

  c(qnew[1:(npop-1),],t(fnew))
}

emUnknown<-function(x,like,npop=2){
  ind<-length(like)
  sites<-ncol(like[[1]])
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  F<-1-t(f)
  expGeno<-function(x){
     fPart<-F%*%q[,x]
     pp<-t(like[[x]])*cbind(fPart^2,2*fPart*(1-fPart),(1-fPart)^2)
     pp<-pp/rowSums(pp)
     pp[,2]+2*pp[,3]
  }
 # genoMat<-matrix(unlist(geno),ncol=length(geno))
  genoMat<-sapply(1:ind,expGeno)


  aNorm<-matrix(NA,nrow=ind,ncol=sites)
  bNorm<-matrix(NA,nrow=ind,ncol=sites)
  for(tal in 1:sites){
    aNorm[,tal]<-t(q)%*%f[,tal]
    bNorm[,tal]<-t(q)%*%(1-f[,tal])
  }

  updateF<-function(j,k){
    ag<-genoMat[j,]*q[k,]*f[k,j]/aNorm[,j]
    bg<-(2-genoMat[j,])*q[k,]*(1-f[k,j])/bNorm[,j]
    fnew<-sum(ag)/(sum(ag)+sum(bg))
    fnew
  }
##  (f[1,1]<-updateF(1,1))
  updateQ<-function(i,k){
    ag<-genoMat[,i]*q[k,i]*f[k,]/aNorm[i,]
    bg<-(2-genoMat[,i])*q[k,i]*(1-f[k,])/bNorm[i,]
    qnew<-sum(bg+ag)/(2*sites)
    qnew
  }

  
  fnew<-f
  qnew<-q
  for(k in 1:npop)
    for(j in 1:sites)
      fnew[k,j]<-updateF(j,k)

  
  for(k in 1:(npop-1))
    for(i in 1:ind)
      qnew[k,i]<-updateQ(i,k)

  
  c(qnew[1:(npop-1),],t(fnew))
}

splitIn<-function(x,mc){
  s<-seq(0,x,by=x/mc)
  s[length(s)]<-x
 ll<-list()
 for(i in 1:mc)
   ll[[i]]<-(s[i]+1):(s[i+1])
ll
}



if(FALSE){
llike<-numeric(ind*sites*3)
dim(llike)<-c(ind,sites*3)
for(i in 1:ind)
  llike[i,]<-like[[i]]

llike<-t(llike)
}
expGenoCPP_input2<-signature(f="numeric",l="numeric",K="integer",q="numeric",N="integer",M="integer")
#updateQ<-function(i,k){
#    ag<-genoMat[,i]*q[k,i]*f[k,]/aNorm[i,]
#    bg<-(2-genoMat[,i])*q[k,i]*(1-f[k,])/bNorm[i,]
#    qnew<-sum(bg+ag)/(2*sites)
#    qnew
#  }

expGenoCPP_code2<-"
double *bNorm=malloc(sizeof(double)*N[0]*M[0]);
double *aNorm=malloc(sizeof(double)*N[0]*M[0]);
double *fnew=malloc(sizeof(double)*N[0]*M[0]);
double *expG=malloc(sizeof(double)*N[0]*M[0]);
double *qnew=malloc(sizeof(double)*N[0]*K[0]);
for(int i=0;i<N[0];i++){
  for(int j=0;j<M[0];j++){
    double fpart=0;
    for(int k=0;k<K[0];k++)
      fpart=fpart+f[j*K[0]+k]*q[i*K[0]+k];
    aNorm[i+j*N[0]]=1.0/fpart;
    fpart=1-fpart;
    bNorm[i+j*N[0]]=1.0/fpart;
    double pp0=fpart*fpart*        l[i*3*M[0]+j*3+0];
    double pp1=2*(1-fpart)*fpart*  l[i*3*M[0]+j*3+1];
    double pp2=(1-fpart)*(1-fpart)*l[i*3*M[0]+j*3+2];
    double sum=pp0+pp1+pp2;
    expG[j+i*M[0]]=(pp1+2*pp2)/sum;
  }

}
 for(int j=0;j<M[0];j++){
    for(int k=0;k<K[0];k++){
      double sumAG=0;
      double sumBG=0;
      for(int i=0;i<N[0];i++){  
        sumAG=sumAG+expG[j+i*M[0]]*f[j*K[0]+k]*q[i*K[0]+k]*aNorm[i+j*N[0]];
        sumBG=sumBG+(2-expG[j+i*M[0]])*q[i*K[0]+k]*(1-f[j*K[0]+k])*bNorm[i+j*N[0]];
      }
      fnew[j*K[0]+k]=sumAG/(sumAG+sumBG);
    }
  }
for(int i=0;i<N[0];i++){  
   for(int k=0;k<K[0];k++){
     double sumAGBG=0;
     for(int j=0;j<M[0];j++){
        sumAGBG=sumAGBG+expG[j+i*M[0]]*f[j*K[0]+k]*q[i*K[0]+k]*aNorm[i+j*N[0]]+
        (2-expG[j+i*M[0]])*q[i*K[0]+k]*(1-f[j*K[0]+k])*bNorm[i+j*N[0]];
     }
     qnew[i*K[0]+k]=sumAGBG/(2*M[0]);
   }
 }

for(int k=0;k<K[0];k++)
  for(int j=0;j<M[0];j++)
    f[j*K[0]+k]=fnew[j*K[0]+k];

for(int k=0;k<K[0];k++)
  for(int i=0;i<N[0];i++)
    q[i*K[0]+k]=qnew[i*K[0]+k];

free(aNorm);
free(fnew);
free(qnew);
free(bNorm);
free(expG);
"
fns2 <- cfunction( list(eG=expGenoCPP_input2), 
                   list(expGenoCPP_code2), 
                   convention=".C", cxxargs="-O3", cppargs="-O3",language="C")
ffun2<-fns2[["eG"]]
#system.time(ff2<-ffun2(f=f,l=llike,K=npop,q=q,ind,sites))

#str(ff2)

####################################3


#loglikeUnKnown3<-function(q,f,llike,npop=2){
#  F<-1-t(f)
#  loglikeIndlike<-function(lPart,qPart){
#    fPart<-F%*%qPart
#    sum(log((fPart^2*lPart[1,]+2*fPart*(1-fPart)*lPart[2,]+(1-fPart)^2*lPart[3,])))
#  }
#  loglike<-0
#  for(tal in 1:ncol(llike))
#    loglike<-loglike+loglikeIndlike(matrix(llike[,tal],3),q[,tal])
#  return(loglike)
#}





getLikeCPP_input3<-signature(f="numeric",l="numeric",K="integer",q="numeric",N="integer",M="integer",logLike="numeric")

getLikeCPP_code3<-"
for(int i=0;i<N[0];i++){
  for(int j=0;j<M[0];j++){
    double fpart=0;
    for(int k=0;k<K[0];k++)
      fpart=fpart+f[j*K[0]+k]*q[i*K[0]+k];


    double pp0=fpart*fpart*        l[i*3*M[0]+j*3+0];
    double pp1=2*(1-fpart)*fpart*  l[i*3*M[0]+j*3+1];
    double pp2=(1-fpart)*(1-fpart)*l[i*3*M[0]+j*3+2];
    double sum=pp0+pp1+pp2;
    logLike[0]=logLike[0]+log(sum);
  }

}
"
fns3 <- cfunction( list(logLike=getLikeCPP_input3), 
                   list(getLikeCPP_code3), 
                   convention=".C", cxxargs="-O3", cppargs="-O3",language="C")
likeCpp<-fns3[["logLike"]]
#system.time(ff2<-ffun2(f=f,l=llike,K=npop,q=q,ind,sites))

#str(ff2)


########################################
expGenoCPP_input<-signature(f="numeric",l="numeric",K="integer",q="numeric",N="integer",M="integer",res="numeric")

# expGeno<-function(x,like,q,freq){
#     fPart<-freq%*%q[,x]
#     pp<-t(like[[x]])*cbind(fPart^2,2*fPart*(1-fPart),(1-fPart)^2)
#     pp<-pp/rowSums(pp)
#     pp[,2]+2*pp[,3]
#  }
expGenoCPP_code<-"

for(int i=0;i<M[0];i++){

  double fpart=0;
  for(int k=0;k<K[0];k++)
    fpart=fpart+f[i+k*M[0]]*q[k];
  double pp0=fpart*fpart*l[i*3+0];
  double pp1=2*(1-fpart)*fpart*l[i*3+1];
  double pp2=(1-fpart)*(1-fpart)*l[i*3+2];
  double sum=pp0+pp1+pp2;

  res[i]=(pp1+2*pp2)/sum;
}
"
fns <- cfunction( list(eG=expGenoCPP_input), 
                   list(expGenoCPP_code), 
                   convention=".C")
fffun<-fns[["eG"]]

#ff<-ffun(f=F,l=like[[1]],K=npop,q=q[,1],ind,sites,numeric(sites))
#str(ff)


emUnknown3<-function(x,llike,npop=2){
  minErr<-1e-10
  kk<<-x
  ind<-ncol(llike)
  sites<-nrow(llike)/3
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  res<-ffun2(f=f,l=llike,K=npop,q=q,ind,sites)
  dim(res$q)<-c(npop,ind)
  dim(res$f)<-c(npop,sites)
  res$f[res$f<minErr]<-minErr
  res$f[res$f>1-minErr]<-1-minErr
  res$q[res$q<minErr]<-minErr
  res$q[res$q>1-minErr*2]<-1-minErr*2
  res$q<-t(t(res$q)/colSums(res$q))
  
  xnew<-c(res$q[1:(npop-1),],t(res$f))
#  barplot(res$q[1:(npop-1),]-x[1:(ind*(npop-1))],beside=T)
  xnew
}

emUnknown2<-function(x,like,npop=2){
  #x<-c(qTrue[2,],cbind(f1,f2))
  kk<<-x
  ind<-length(like)
  sites<-ncol(like[[1]])
  q<-matrix(x[1:(ind*(npop-1))],nrow=npop-1)
  q<-rbind(q,1-apply(q,2,sum))
  f<-matrix(x[-(1:(ind*(npop-1)))],nrow=npop,by=T)
  #F<-1-t(f)
  ffout<-ffun2(f=f,l=llike,K=npop,q=q,ind,sites,numeric(sites*ind),numeric(sites*ind),numeric(sites*ind))
  genoMat<-ffout$res
  dim(genoMat)<-c(sites,ind)
  aNorm<-ffout$aNorm
  bNorm<-ffout$bNorm
  dim(aNorm)<-c(ind,sites)
  dim(bNorm)<-c(ind,sites)
  #expGeno<-function(x,like,q,freq){
  #   fPart<-freq%*%q[,x]
  #   pp<-t(like[[x]])*cbind(fPart^2,2*fPart*(1-fPart),(1-fPart)^2)
  #   pp<-pp/rowSums(pp)
  #   pp[,2]+2*pp[,3]
  #}
  
 # genoMat<-sapply(1:ind,expGeno,like=like,q=q,freq=F)
  #genoMat<-unlist(mclapply(splitIn(ind,mc.cores),function(y) sapply(y,expGeno,like=like,q=q,freq=F),mc.cores=mc.cores))
 # dim(genoMat)<-c(sites,ind)
  #genoMat<-sapply(1:ind,function(x) ffun(f=F,l=like[[x]],K=npop,q=q[,x],ind,sites,numeric(sites))$res)
 # aNorm<-sapply(1:sites,function(tal) t(q)%*%f[,tal])
 # bNorm<-sapply(1:sites,function(tal) t(q)%*%(1-f[,tal]))
  #aNorm<-matrix(unlist(mclapply(1:sites,function(tal) t(q)%*%f[,tal],mc.cores=mc.cores)),nrow=ind)
  #bNorm<-matrix(unlist(mclapply(1:sites,function(tal) t(q)%*%(1-f[,tal]),mc.cores=mc.cores)),nrow=ind)

 # aa<-mclapply(splitIn(sites,mc.cores),function(y) sapply(y,function(tal) t(q)%*%f[,tal]),mc.cores=mc.cores)
 # aNorm<-matrix(unlist(aa),nrow=ind)

 # bb<-mclapply(splitIn(sites,mc.cores),function(y) sapply(y,function(tal) t(q)%*%(1-f[,tal])),mc.cores=mc.cores)
 # bNorm<-matrix(unlist(bb),nrow=ind)

  
  updateF<-function(j,k){
    ag<-genoMat[j,]*q[k,]*f[k,j]/aNorm[,j]
    bg<-(2-genoMat[j,])*q[k,]*(1-f[k,j])/bNorm[,j]
    fnew<-sum(ag)/(sum(ag)+sum(bg))
    fnew
  }
##  (f[1,1]<-updateF(1,1))
  updateQ<-function(i,k){
    ag<-genoMat[,i]*q[k,i]*f[k,]/aNorm[i,]
    bg<-(2-genoMat[,i])*q[k,i]*(1-f[k,])/bNorm[i,]
    qnew<-sum(bg+ag)/(2*sites)
    qnew
  }

  
  fnew<-f
  qnew<-q
  for(k in 1:npop)
    fnew[k,]<-unlist(mclapply(1:sites,updateF,k=k,mc.cores=mc.cores))

 
  
  for(k in 1:(npop-1))
    qnew[k,]<-unlist(mclapply(1:ind,updateQ,k=k,mc.cores=mc.cores))


  qqNew=qnew[1:(npop-1),]
  if(any(is.na(qqNew))){
    cat("QQ old\n")
    print(q)
    print(f)
    cat("QQ new\n")
      print(qqNew)
  }
  fnew[fnew<0.0001]<-0.0001
   c(qnew[1:(npop-1),],t(fnew))
}
#system.time(h<-emUnknown2(x,like=like,npop=npop))



#system.time(replicate(1,res<-emUnknown2(c(qTrue[2,],cbind(f1,f2)),like)))

start<-function(geno,npop=2){
  ind<-length(geno)
  sites<-length(geno[[1]])
  q<-rep(1/npop,(npop-1)*ind)
  f<-apply(matrix(unlist(geno),ncol=length(geno)),1,mean)/2
  c(q,f,f)
}
start2<-function(geno,npop=2){
  ind<-length(geno)
  if(class(geno[[1]])=="numeric")
    sites<-length(geno[[1]])
  else
    sites<-ncol(geno[[1]])
  q<-matrix(runif(npop*ind),ncol=npop)
  q<-t(q/rowSums(q))
  f<-runif(sites*npop)
  c(q[-npop,],f)
}

sq<-function (par, fixptfn, objfn,maxiter=1500,tol=1e-7,step.max0=1,mstep=4,objfn.inc=1,llike,npop,alpha.max=Inf)
#   maxiter=1500;tol=1e-7;step.min=1;step.max0=1;mstep=4;objfn.inc=1
#  par<-start;fixptfn=emUnknown3;objfn=loglikeUnKnownWrap3;llike=llike;tol=0.0001;npop=npop
{
  step.max=step.max0;    
  p <- par
  
  lold <- objfn(p,llike=llike,npop=npop)
  for(iter in 1:maxiter) {
    p1 <- fixptfn(p,llike=llike,npop=npop)
    q1 <- p1 - p
    sr2 <- sum(q1^2)

    if (sqrt(sr2) < tol) 
      break

    p2 <- fixptfn(p1,llike=llike,npop=npop)
    q2 <- p2 - p1
    sq2 <- sqrt(sum(q2^2))
    if (sq2 < tol) 
      break
    sv2 <- sum((q2 - q1)^2)
    srv <-  sum(q1*(q2-q1))
    alpha <- sqrt(sr2/sv2)
    alpha <- max(1, min(step.max, alpha))
    p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1)

   # if(any(p.new<0 | p.new>1)){
   #   print(    range(p.new))
   #   kk<<-p.new
   #   stop()
   # }
    if (abs(alpha - 1) > 0.01) {
      p.new <- fixptfn(p.new,llike=llike,npop=npop)
     
    }
    lnew <- objfn(p.new,llike=llike,npop=npop)
    #lp2 <- objfn(p2,llike=llike,npop=npop)
    if ( (lnew > lold + objfn.inc)) {
      p.new <- p2
      lnew <- objfn(p2,llike=llike,npop=npop)
      if ((alpha - step.max)<1e6) 
        step.max <- max(step.max0, step.max/mstep)
      alpha <- 1
      cat("no ")
    }
    cat(alpha,"\t",step.max,"\t",lnew,"\n")
    if (alpha == step.max) 
      step.max <- mstep * step.max
    # barplot(p[1:40]-p.new[1:40],main="iter")
    p <- p.new
    lold<-lnew
  }

  return(list(par = p, value.objfn = lnew, iter = iter))
}
