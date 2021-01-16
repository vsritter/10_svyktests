compCIF <- function(time, cens,
                    group = rep(1, length(time)),
                    weights = NULL, return_si = F){
  # this function compares the CIF-s of two groups      #
  # it is based on the test presented in the paper      #
  # from Statistics in Medicine, 1993, by Pepe and Mori #
  ttau=unlist(by(time,group,max))
  tau=min(ttau)
  
  if(is.null(weights)) {
    ng=table(group)
  } else {
    ng=tapply(weights, group, sum)
  }
  lg=labels(ng)$group
  
  if(is.null(weights)) {
    dd=table(time,cens,group)
  } else {
    dd=table(time,cens,group)*weights[order(time)]
  }
  
  if (dim(dd)[3]!=2)
    stop("Pepe-Mori test is for two groups")
  
  if (dim(dd)[2]>3)
    stop("All competing risks should be grouped under code 2")
  
  if (dim(dd)[2]==2){
    ddd=array(0,dim=c(dim(dd)[1],3,2))
    
    if (sum(cens==0)==0){
      ddd[,2:3,]=dd
    }
    
    if (sum(cens==1)==0){
      warning('There are only competing risks')
      ddd[,1,]=dd[,1,]
      ddd[,3,]=dd[,2,]
    }
    
    if (sum(cens==2)==0){
      warning('There are no competing risks')
      ddd[,1:2,]=dd
    }
    
    dd=ddd
  }
  
  if (dim(dd)[2]<2)
    stop('Either all observations are censored or\n
      there is only one type of event and no censor
      observations')
  
  tt=sort(unique(time))
  tt1=c(tt[2:length(tt)],NA)
  
  nt=table(tt1<=tau)[2]
  deltat=tt1-tt
  dd=dd[1:nt,,]
  deltat=deltat[1:nt]
  tt=tt[1:nt]
  
  for (i in 1:2)
  {
    dd1=dd[,,i] ## table of censor, ev, cr
    dd2=apply(dd1,1,sum) ## sum for each time point
    nrisk=ng[i]-cumsum(c(0,dd2[1:(nt-1)]))
    dev=dd1[,2]
    dcr=dd1[,3]
    dcens=dd1[,1]+dcr
    dall=dev+dcr
    si=(nrisk-dall)/nrisk
    s=cumprod(si)
    sminus=c(1,s[1:(length(s)-1)])
    fi=dev/nrisk*sminus
    f=cumsum(fi)
    fcri=dcr/nrisk*sminus
    fcr=cumsum(fcri)
    Ci=(nrisk-dcens)/nrisk
    C=cumprod(Ci)
    Cminus=c(1,C[1:(length(C)-1)])
    
    h = nrisk/sminus
    
    if (i==1) {
      nrisk1=nrisk
      s1=s
      f1=f
      fcr1=fcr
      C1=C
      C1minus=Cminus
      dall1=dall
      dev1=dev
      
      h1=h
    }
    
    if(i==2){
      nrisk2=nrisk
      s2=s
      f2=f
      fcr2=fcr
      C2=C
      C2minus=Cminus
      dall2=dall
      dev2=dev
      
      h2=h
    }
  }
  
  wi=C1minus*C2minus*sum(ng)/(ng[1]*C1minus+ng[2]*C2minus)
  si=wi*(f1-f2)*deltat
  s=sqrt(ng[1]*ng[2]/sum(ng))*sum(si)

  phii = sqrt(ng[1]*ng[2]/sum(ng))*wi*(dev1/h1-dev2/h2)
  
  # for group 1
  temp=wi*(1-f1)*deltat
  temp[is.na(temp)]=0
  t1parti=rev(cumsum(rev(temp)))
  temp=wi*f1*deltat
  temp[is.na(temp)]=0
  t2i=rev(cumsum(rev(temp)))
  temp=wi*deltat
  temp[is.na(temp)]=0
  t3i=fcr1*rev(cumsum(rev(temp)))
  t1i=t1parti-t3i
  sigma1i=(dev1*t1i^2+(dall1-dev1)*t2i^2)/(nrisk1*(nrisk1-1))
  sigma1=sum(sigma1i,na.rm=T)
  
  # for group 2
  temp=wi*(1-f2)*deltat
  temp[is.na(temp)]=0
  t1parti=rev(cumsum(rev(temp)))
  temp=wi*f2*deltat
  temp[is.na(temp)]=0
  t2i=rev(cumsum(rev(temp)))
  temp=wi*deltat
  temp[is.na(temp)]=0
  t3i=fcr2*rev(cumsum(rev(temp)))
  t1i=t1parti-t3i
  sigma2i=(dev2*t1i^2+(dall2-dev2)*t2i^2)/(nrisk2*(nrisk2-1))
  sigma2=sum(sigma2i,na.rm=T)
  
  sigma=ng[1]*ng[2]*(sigma1+sigma2)/sum(ng)
  z=s^2/sigma
  pvalue=1-pchisq(z,1)
  # r=data.frame(zstat = s, var = sigma,
  #              chisq = z, p = pvalue)
  r=data.frame(s = s, v = sigma)
  row.names(r)=''
  
  if(return_si){
    return(list(
      as.vector(phii),
      ng[1]*sum(ng)/ng[2],
      s, sigma))
  } else {
    return(r)
  }
  # return(data.frame(sminus1, sminus2, dev1, dev2, nrisk1, nrisk2))
}



