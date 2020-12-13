`mfma` <-
function(y,k1,k2,r,it=15,eps=0.001,init=NULL,scaling=FALSE,seed=7,model.names=c(1,1))
{
ptm <- proc.time()
set.seed(seed)

if (scaling) y<-scale(y)
numobs<-nrow(y)
p<-ncol(y)
ybar<- apply(y, 2, mean)
y<-scale(y, ybar, scale=FALSE) 
lik<--100000000000
muf=matrix(0,p,k1)
### operazioni preliminari ed inizializzazione dei parametri
### l'inizializzazione è casuale se non altrimenti specificata per mu e sigma
### per i pesi l'inizializzazione viene fatta sulla base della partizione k medie

#output<-hc(modelName = "VVV", data = y)
#if (k1>1) s1<-hclass(output,k1) else s1<-rep(1,numobs)

if (k1>1) s1=kmeans(y,k1,iter.max = 1000)$cl else s1=rep(1,numobs)
for  (i in 1:k1) if ((table(s1)[i])<2) s1[sample(1:numobs,2,replace=FALSE)]=i

psi<-array(0,c(k1,p,p))
H<-array(0,c(k1,p,r))


if (is.null(init)) {for (j in 1:k1) { stima=try(factanal(y[s1==j,],r,rotation="none"),silent=TRUE) 
                                          if (is.character(stima)) {psi[j,,]<-0.1*diag(p)
                                           H[j,,]<-matrix(runif(p*r),p,r)}

                                          if (!is.character(stima)) {psi[j,,]<-diag(stima$uniq)
                                                                         H[j,,]<-stima$load}

                            }} else {psi<-init$psi
                                    H=init$H}
                                               
                                               
if (is.null(init$muf)) for (j in 1:k1) muf[,j]=colMeans(y[s1==j,])
if (is.null(init$w1)) w1<-matrix(table(s1)/numobs) else w1<-init$w1

##############################################################################

z=matrix(factanal(y,r,scores="Bartlett")$scores,numobs)
if (k2>1) s2=kmeans(z,k2)$cl else s2=rep(1,numobs)
for  (i in 1:k2) if ((table(s2)[i])<2) s2[sample(1:numobs,2,replace=FALSE)]=i

mu=matrix(0,k2,r)                     
sigma<-array(0,c(k2,r,r))

for (i in 1:k2) {
                mu[i,]<-colMeans(matrix(z[s2==i,],ncol=r))
                sigma[i,,]<-var(z[s2==i,])
                }

if (is.null(init$w2)) w2<-matrix(table(s2)/numobs) else w2<-init$w2
#mu[k2,]<--t(w2[1:k2-1])%*%mu[1:k2-1,]/w2[k2]

w1<-matrix(w1)
w2<-matrix(w2)

out<-mfma.em.alg(y,numobs,r,k1,k2,p,it,H,w1,w2,mu,sigma,eps,psi,lik,muf,model.names)

likelihood<-out$likelihood
ps1.y<-out$ps1.y
ps2.y<-out$ps2.y
sigma<-out$sigma
H<-out$H
w2<-out$w2
w1<-out$w1
mu<-out$mu
psi<-out$psi
ps1s2.y=out$ps1s2.y
s12.bis=out$s12
likelihood<-matrix(likelihood[!likelihood==0])
s1=out$s1
s2=out$s2
muf=out$muf
s12=s2+((s1-1)*s1+(s1-1)*1)


#z=matrix(0,numobs,r)
#for (j in 1:k1) z[s1==j,]<-t(ginv(t(H[j,,])%*%ginv(psi[j,,])%*%H[j,,])%*%t(H[j,,])%*%ginv(psi[j,,])%*%t(matrix(y[s1==j,],ncol=p)))

if (model.names[1]==1) h1=(k1-1)+(p*r)*k1+p*k1-(k1*r*(r-1)/2)+(k1-1)*p
if (model.names[2]==1) h2=(k2-1)+r*(k2-1)+(r*(r+1)/2)*(k2-1)

if (model.names[1]==2) h1=(k1-1)+(p*r)*k1+p-(k1*r*(r-1)/2)+(k1-1)*p
if (model.names[2]==2) h2=(k2-1)+r*(k2-1)

if (model.names[1]==3) h1=(k1-1)+(p*r)*k1+k1-(k1*r*(r-1)/2)+(k1-1)*p
if (model.names[2]==3) h2=(k2-1)+r*(k2-1)+(k2-1)



h=h1+h2
lik<-likelihood[length(likelihood)]
bic<--2*lik+h*log(numobs)
aic<--2*lik+2*h
EN<-entr(matrix(ps1s2.y,k2*k1,numobs))
clc<--2*lik+2*EN
icl.bic<--2*lik+2*EN+h*log(numobs)

output<-list(h=h,k1=k1,H=H,lik=likelihood,w2=w2,w1=w1,mu=mu,
             psi=psi,k2=k2,sigma=sigma,p=p,numobs=numobs,ph.y=ps1s2.y,
             scaling=scaling,r=r,s2=s2,s1=s1,bic=bic,aic=aic,s12=s12.bis,
             clc=clc,icl.bic=icl.bic,cptime=proc.time()-ptm,seme=seed,muf=muf,h=h)
invisible(output)
}

