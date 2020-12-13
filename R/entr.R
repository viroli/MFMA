`entr` <-
function(z)
{
numobs<-nrow(z)
numg<-ncol(z)
temp<-0
z<-ifelse(z==0,z+0.000000000000000000000001,z)
for (i in 1:numg) for (j in 1:numobs) temp<-temp+(z[j,i]*log(z[j,i]))
return(-temp)
}

