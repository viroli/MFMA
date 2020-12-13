`mfa` <-
function(y,k,r,it=15,eps=0.001,init=NULL,scaling=FALSE,seed=7,model.names=1)
{
ptm <- proc.time()
set.seed(seed)

k2=1
model.names=c(model.names,1)
if (scaling) y<-scale(y)
numobs<-nrow(y)
p<-ncol(y)
ybar<- apply(y, 2, mean)
y<-scale(y, ybar, scale=FALSE)
lik<--100000000000
muf=matrix(0,p,k)
### operazioni preliminari ed inizializzazione dei parametri
### l'inizializzazione ? casuale se non altrimenti specificata per mu e sigma
### per i pesi l'inizializzazione viene fatta sulla base della partizione k medie

#output<-hc(modelName = "VVV", data = y)
#if (k>1) s1<-hclass(output,k) else s1<-rep(1,numobs)

if (k>1) s1=kmeans(y,k,iter.max = 1000)$cl else s1=rep(1,numobs)
for  (i in 1:k) if ((table(s1)[i])<2) s1[sample(1:numobs,2,replace=FALSE)]=i

psi<-array(0,c(k,p,p))
H<-array(0,c(k,p,r))


if (is.null(init)) {for (j in 1:k) { stima=try(factanal(y[s1==j,],r,rotation="none"),silent=TRUE)
                                          if (is.character(stima)) {psi[j,,]<-0.1*diag(p)
                                           H[j,,]<-matrix(runif(p*r),p,r)}

                                          if (!is.character(stima)) {psi[j,,]<-diag(stima$uniq)
                                                                         H[j,,]<-stima$load}

                            }} else {psi<-init$psi
                                    H=init$H}


if (is.null(init$muf)) for (j in 1:k) muf[,j]=colMeans(y[s1==j,])
if (is.null(init$w1)) w1<-matrix(table(s1)/numobs) else w1<-init$w1

##############################################################################

z=matrix(factanal(y,r,scores="Bartlett")$scores,numobs)
s2=rep(1,numobs)
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

out<-mfma.em.alg(y,numobs,r,k,k2,p,it,H,w1,w2,mu,sigma,eps,psi,lik,muf,model.names)

likelihood<-out$likelihood
ps1.y<-out$ps1.y
H<-out$H
w1<-out$w1
psi<-out$psi
likelihood<-matrix(likelihood[!likelihood==0])
s1=out$s1
muf=out$muf
ps1s2.y=out$ps1s2.y

#z=matrix(0,numobs,r)
#for (j in 1:k) z[s1==j,]<-t(ginv(t(H[j,,])%*%ginv(psi[j,,])%*%H[j,,])%*%t(H[j,,])%*%ginv(psi[j,,])%*%t(matrix(y[s1==j,],ncol=p)))

if (model.names[1]==1) h1=(k-1)+(p*r)*k+p*k-(k*r*(r-1)/2)+(k-1)*p
if (model.names[1]==2) h1=(k-1)+(p*r)*k+p-(k*r*(r-1)/2)+(k-1)*p
if (model.names[1]==3) h1=(k-1)+(p*r)*k+k-(k*r*(r-1)/2)+(k-1)*p
h2=(k2-1)+r*(k2-1)+(r*(r+1)/2)*(k2-1)



h=h1+h2
lik<-likelihood[length(likelihood)]
bic<--2*lik+h*log(numobs)
aic<--2*lik+2*h
EN<-entr(matrix(ps1s2.y,k2*k,numobs))
clc<--2*lik+2*EN
icl.bic<--2*lik+2*EN+h*log(numobs)

output<-list(h=h,k=k,H=H,lik=likelihood,w=w1,psi=psi,p=p,numobs=numobs,
             scaling=scaling,r=r,cl=s1,bic=bic,aic=aic,clc=clc,icl.bic=icl.bic,
             cptime=proc.time()-ptm,seme=seed,mu=muf,h=h,ps.y=ps1.y)
invisible(output)
}

