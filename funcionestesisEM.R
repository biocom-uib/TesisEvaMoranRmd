#### Tests

#Test para diferencias de RA

RA.RA=function(M,conf.level=0.95){
  beta=c(1,-1,-1,1)
  B=sum(beta)
  e=M[1,]
  n=colSums(M)
  p=e/n
  q=1-2*p
  N=sum(n)
  L.b=sum(beta*p)
  T=L.b+n[2]+n[3]
  T0=4*L.b^2/sum(beta^2/n)
  if (L.b>0){T1=T*L.b/2}
  if (L.b<0){T1=(N-T)*L.b/(-2)}
  if (L.b==0){T1=0}
  y=function(x){N-sum(sqrt(n^2+((2*n*beta*q)/L.b)*x+(beta^2/L.b^2)*x^2))}
  z2=rootSolve::multiroot(y,start=(T0+T1)/2,maxiter=1000,rtol = 1e-8, atol = 1e-10, ctol = 1e-10)$root
  p.val=1-pchisq(z2,1)
  zz=qnorm((1-conf.level)/2)^2
  y=function(x){N+2*zz -2*L.b*x-sum(sqrt(n^2+2*n*beta*q*x+x^2))}
  C1=rootSolve::multiroot(y,start=-100,maxiter=1000,rtol = 1e-8, atol = 1e-10, ctol = 1e-10)$root
  C2=rootSolve::multiroot(y,start=100,maxiter=1000,rtol = 1e-8, atol = 1e-10, ctol = 1e-10)$root
  IC=L.b-zz/c(C2,C1)
  Res=rbind(c(p.val,IC))
  colnames(Res)=c("p-valor","IC-inferior","IC-superior")
  return(Res)
}


#Test para cociente de RR

RR.RR=function(M,conf.level=0.95){
  e=M[1,]
  n=colSums(M)
  p=e/n
  
  RR1=p[1]/p[2]
  RR2=p[3]/p[4]
  se.lnRR1=sqrt(1/e[1]+1/e[2]-1/n[1]-1/n[2])
  se.lnRR2=sqrt(1/e[3]+1/e[4]-1/n[3]-1/n[4])
  
  RRR=RR1/RR2
  se.lnRRR=sqrt(se.lnRR1^2+se.lnRR2^2)
  
  IC.1=exp(log(RRR)-se.lnRRR*qnorm((1+conf.level)/2))
  IC.2=exp(log(RRR)+se.lnRRR*qnorm((1+conf.level)/2))
  p.val=2-2*pnorm(abs(log(RRR))/se.lnRRR)
  Res=rbind(c(RRR,p.val,IC.1,IC.2))
  colnames(Res)=c("RRR", "p-valor","IC-inferior","IC-superior")
  return(Res)
}

# Test para comparar OR

OR.OR=function(M,conf.level=0.95){
  M1=M[,1:2]
  M2=M[,3:4]
  LOR1=log(fisher.test(M1)$estimate)
  LOR2=log(fisher.test(M2)$estimate)
  Z=abs(LOR1-LOR2)/sqrt(sum(1/M))
  p.val=2*(1-pnorm(Z))
  IC=exp(LOR1-LOR2+qnorm((1+conf.level)/2)*sqrt(sum(1/M))*c(-1,1))
  Res=rbind(c(p.val,IC))
  colnames(Res)=c("p-valor","IC-inferior","IC-superior")
  return(Res)
}

# Test de tendencia

Prop.trend.test=function(x,X){
  chisq.test(rbind(x,X))$p.value  
}


#### Tablas y gráficos


# Tablas 2d Descripción casos-controles
Tabla.DMG=function(I,NI,L1,L2,r=4){
  n_I=length(I)  
  n_NI=length(NI)  
  NA_I=length(I[is.na(I)])
  NA_NI=length(NI[is.na(NI)])
  
  EE=rbind(as.vector(table(I))[2:1],
           as.vector(table(NI))[2:1])
  FT=fisher.test(EE)
  OR=round(FT$estimate,2)
  IC=round(FT$conf.int,2)
  
  EEExt=rbind(c(as.vector(table(I))[2:1],NA_I),
              c(round(100*as.vector(table(I))[2:1]/(n_I-NA_I),1),NA),
              c(as.vector(table(NI))[2:1],NA_NI),
              c(round(100*as.vector(table(NI))[2:1]/(n_NI-NA_NI),1),NA),
              c(round(FT$estimate,2),NA,NA),
              c(round(FT$conf.int[1],2),NA,NA),
              c(round(FT$conf.int[2],2),NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Infectadas (N)", "Infectadas (%)","No infectadas (N)", "No infectadas (%)","OR","Extr. Inf. IC","Extr. Sup. IC","p-valor")
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}
#######
#######

Tabla.DMGm=function(I,NI,L,r=4){
  l=length(L)
  n_I=length(I)  
  n_NI=length(NI)  
  NA_I=length(I[is.na(I)])
  NA_NI=length(NI[is.na(NI)])
  EE=rbind(as.vector(table(I)),
           as.vector(table(NI)))
  
  pp=rep(0,dim(EE)[2])
  for (i in 1:dim(EE)[2]){
    pp[i]=prop.test(c(EE[1,i],sum(EE[1,-i])),
                    c(sum(EE[,i]),sum(EE[,-i])))$p.value
  }
  pp=round(p.adjust(pp,method="bonferroni"),r)
  OR=rep(0,dim(EE)[2])
  CI1=rep(0,dim(EE)[2])
  CI2=rep(0,dim(EE)[2])
  
  for (i in 1:dim(EE)[2]){
    FT=fisher.test(rbind(c(EE[1,i],sum(EE[1,-i])),
                         c(EE[2,i],sum(EE[2,-i]))))
    OR[i]=round(FT$estimate,2)
    CI1[i]=round(FT$conf.int[1],2)
    CI2[i]=round(FT$conf.int[2],2)}
  
  
  EEExt=rbind(c(as.vector(table(I)),NA_I),
              c(round(100*as.vector(table(I))/(n_I-NA_I),1),NA),
              c(as.vector(table(NI)),NA_NI),
              c(round(100*as.vector(table(NI))/(n_NI-NA_NI),1),NA),
              c(OR,NA),
              c(CI1,NA),
              c(CI2,NA),
              c(pp,NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Infectadas (N)", "Infectadas (%)","No infectadas (N)", "No infectadas (%)","OR","Extr. Inf. IC","Extr. Sup. IC", "p-valor ajustado")
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}

#########
#########

Tabla.DMGC=function(I,NI,L1,L2,r=4){
  n_I=length(I)  
  n_NI=length(NI)  
  NA_I=length(I[is.na(I)])
  NA_NI=length(NI[is.na(NI)])
  
  EE=rbind(as.vector(table(I))[2:1],
           as.vector(table(NI))[2:1])
  FT=fisher.test(EE)
  PT=prop.test(EE[,1],rowSums(EE))
  
  RA=round(PT$estimate[1]-PT$estimate[2],4)
  ICRA=round(PT$conf.int,4)
  RR=round(PT$estimate[1]/PT$estimate[2],2)
  ICRR=round(RelRisk(EE,conf.level=0.95,method="score")[2:3],2)
  
  EEExt=rbind(c(as.vector(table(I))[2:1],NA_I),
              c(round(100*as.vector(table(I))[2:1]/(n_I-NA_I),1),NA),
              c(as.vector(table(NI))[2:1],NA_NI),
              c(round(100*as.vector(table(NI))[2:1]/(n_NI-NA_NI),1),NA),
              c(RA,NA,NA),
              c(ICRA[1],NA,NA),
              c(ICRA[2],NA,NA),
              c(RR,NA,NA),
              c(ICRR[1],NA,NA),
              c(ICRR[2],NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Infectadas (N)", "Infectadas (%)","No infectadas (N)", "No infectadas (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}
####
####
Tabla.DMGCr=function(I,NI,L1,L2,ni,nni,r=4){
  NA_I=length(I[is.na(I)])
  NA_NI=length(NI[is.na(NI)])
  
  EE=rbind(as.vector(table(I))[2:1],
           as.vector(table(NI))[2:1])
  FT=fisher.test(EE)
  PT=prop.test(EE[,1],rowSums(EE))
  
  RA=round(PT$estimate[1]-PT$estimate[2],4)
  ICRA=round(PT$conf.int,4)
  RR=round(PT$estimate[1]/PT$estimate[2],2)
  ICRR=round(RelRisk(EE,conf.level=0.95,method="score")[2:3],2)
  
  EEExt=rbind(c(as.vector(table(I))[2:1],NA_I),
              c(round(100*as.vector(table(I))[2:1]/(ni-NA_I),1),NA),
              c(as.vector(table(NI))[2:1],NA_NI),
              c(round(100*as.vector(table(NI))[2:1]/(nni-NA_NI),1),NA),
              c(RA,NA,NA),
              c(ICRA[1],NA,NA),
              c(ICRA[2],NA,NA),
              c(RR,NA,NA),
              c(ICRR[1],NA,NA),
              c(ICRR[2],NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Infectadas (N)", "Infectadas (%)","No infectadas (N)", "No infectadas (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}

######
######
Tabla.DMGCm=function(I,NI,L,r=4){
  l=length(L)
  n_I=length(I)  
  n_NI=length(NI)  
  NA_I=length(I[is.na(I)])
  NA_NI=length(NI[is.na(NI)])
  
  EE=rbind(as.vector(table(I)),
           as.vector(table(NI)))
  pp=rep(0,dim(EE)[2])
  RA=rep(0,dim(EE)[2])
  RR=rep(0,dim(EE)[2])
  CIRA1=rep(0,dim(EE)[2])
  CIRA2=rep(0,dim(EE)[2])
  CIRR1=rep(0,dim(EE)[2])
  CIRR2=rep(0,dim(EE)[2])
  
  for (i in 1:dim(EE)[2]){
    EE.temp=rbind(c(EE[1,i],sum(EE[1,-i])),
                  c(EE[2,i],sum(EE[2,-i])))
    FT=fisher.test(EE.temp)
    PT=prop.test(EE.temp[,1],rowSums(EE.temp))
    pp[i]=FT$p.value
    RA[i]=round(PT$estimate[1]-PT$estimate[2],4)
    RR[i]=round(PT$estimate[1]/PT$estimate[2],3)
    CIRA1[i]=round(PT$conf.int,4)[1]
    CIRA2[i]=round(PT$conf.int,4)[2]
    CIRR1[i]=round(RelRisk(EE.temp,conf.level=0.95,method="wald"),3)[2]
    CIRR2[i]=round(RelRisk(EE.temp,conf.level=0.95,method="wald"),3)[3]
  }
  pp=p.adjust(pp,method="bonferroni")
  
  EEExt=rbind(c(as.vector(table(I)),NA_I),
              c(round(100*as.vector(table(I))/(n_I-NA_I),1),NA),
              c(as.vector(table(NI)),NA_NI),
              c(round(100*as.vector(table(NI))/(n_NI-NA_NI),1),NA),
              c(RA,NA),
              c(CIRA1,NA),
              c(CIRA2,NA),
              c(RR,NA),
              c(CIRR1,NA),
              c(CIRR2,NA),
              c(round(pp,r),NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Infectadas (N)", "Infectadas (%)","No infectadas (N)", "No infectadas (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}
####
####
# Tablas 2d Descripción sintomatología

Tabla.DMCasos=function(I,L1,L2,r=6){
  I=factor(I)
  IA=I[Casos$SINTOMAS_DIAGNOSTICO=="Asintomática"]
  IL=I[Casos$SINTOMAS_DIAGNOSTICO=="Leve"]
  IG=I[Casos$SINTOMAS_DIAGNOSTICO=="Grave"]
  NA_A=length(IA[is.na(IA)])
  NA_L=length(IL[is.na(IL)])
  NA_G=length(IG[is.na(IG)])
  
  EEExt=cbind(c(as.vector(table(IA))[2:1],NA_A),
              c(round(100*as.vector(prop.table(table(IA))[2:1]),1),NA),
              c(as.vector(table(IL))[2:1],NA_L),
              c(round(100*as.vector(prop.table(table(IL))[2:1]),1),NA),
              c(as.vector(table(IG))[2:1],NA_G),
              c(round(100*as.vector(prop.table(table(IG))[2:1]),1),NA)
  )
  
  test=chisq.test(EEExt[1:2,c(1,3,5)])
  if(any(test$expected< 5))
  {
    p.val=chisq.test(EEExt[1:2,c(1,3,5)],simulate.p.value=TRUE,B=10000)$p.value
    tipus="Montecarlo"} else
    {
      p.val=test$p.value
      tipus="Paramétrico"
    }
  p.val=round(p.val,r)
  EEExt=cbind(data.frame(EEExt),c(p.val,NA,NA),c(tipus,NA,NA))
  rownames(EEExt)=c(L1,L2,"Datos perdidos")
  colnames(EEExt)=c(Columnes.Sint,"p-valor","Tipo")
  
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}
######
#####
Tabla.DMCasosr=function(IA,IL,IG,L1,L2,r=6){
  NA_A=length(IA[is.na(IA)])
  NA_L=length(IL[is.na(IL)])
  NA_G=length(IG[is.na(IG)])
  
  EEExt=cbind(c(as.vector(table(IA))[2:1],NA_A),
              c(round(100*as.vector(prop.table(table(IA))[2:1]),1),NA),
              c(as.vector(table(IL))[2:1],NA_L),
              c(round(100*as.vector(prop.table(table(IL))[2:1]),1),NA),
              c(as.vector(table(IG))[2:1],NA_G),
              c(round(100*as.vector(prop.table(table(IG))[2:1]),1),NA)
  )
  
  test=chisq.test(EEExt[1:2,c(1,3,5)])
  if(any(test$expected< 5))
  {
    p.val=chisq.test(EEExt[1:2,c(1,3,5)],simulate.p.value=TRUE,B=10000)$p.value
    tipus="Montecarlo"} else
    {
      p.val=test$p.value
      tipus="Paramétrico"
    }
  p.val=round(p.val,r)
  EEExt=cbind(data.frame(EEExt),c(p.val,NA,NA),c(tipus,NA,NA))
  rownames(EEExt)=c(L1,L2,"Datos perdidos")
  colnames(EEExt)=c(Columnes.Sint,"p-valor","Tipo")
  
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}

######
#####
Tabla.DMCasosm=function(I,L,r=6){
  l=length(L)
  IA=I[Casos$SINTOMAS_DIAGNOSTICO=="Asintomática"]
  IL=I[Casos$SINTOMAS_DIAGNOSTICO=="Leve"]
  IG=I[Casos$SINTOMAS_DIAGNOSTICO=="Grave"]
  NA_A=length(IA[is.na(IA)])
  NA_L=length(IL[is.na(IL)])
  NA_G=length(IG[is.na(IG)])
  
  EEExt=cbind(c(as.vector(table(IA)),NA_A),
              c(round(100*as.vector(prop.table(table(IA))),1),NA),
              c(as.vector(table(IL)),NA_L),
              c(round(100*as.vector(prop.table(table(IL))),1),NA),
              c(as.vector(table(IG)),NA_G),
              c(round(100*as.vector(prop.table(table(IG))),1),NA)
  )
  p.val=rep(NA,l+1)
  tipus=rep(NA,l+1)
  taula=EEExt[1:l,c(1,3,5)]
  for (i in 1:l){
    X=rbind(taula[i,],colSums(taula)-taula[i,])
    test=chisq.test(X)
    if(any(test$expected< 5))
    {
      p.val[i]=chisq.test(X,simulate.p.value=TRUE,B=10000)$p.value
      tipus[i]="Montecarlo"}
    else
    {
      p.val[i]=test$p.value
      tipus[i]="Paramétrico"
    }
  }
  p.val=round(p.adjust(p.val,method="bonferroni"),r)
  EEExt=cbind(data.frame(EEExt),p.val,tipus)
  rownames(EEExt)=c(L,"Datos perdidos")
  colnames(EEExt)=c(Columnes.Sint,"p-valor","Tipo")
  
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")
}
######
#####
#
# Tabla de comparación Anteparto-Periparto
#
TablaAntePeriOR=function(I,L1,L2,r=6){
  I.AP=I[Casos$PreP==1]
  I.PP=I[Casos$PreP==0]
  n_IAP=length(I.AP)
  n_IPP=length(I.PP)
  NA_IAP=length(I.AP[is.na(I.AP)])
  NA_IPP=length(I.PP[is.na(I.PP)])
  EE=rbind(as.vector(table(I.AP))[2:1],
           as.vector(table(I.PP))[2:1])
  FT=fisher.test(EE)
  OR=round(FT$estimate,2)
  IC=round(FT$conf.int,2)
  
  EEExt=rbind(c(as.vector(table(I.AP))[2:1],NA_IAP),
              c(round(100*as.vector(prop.table(table(I.AP)))[2:1],1),NA),
              c(as.vector(table(I.PP))[2:1],NA_IPP),
              c(round(100*as.vector(prop.table(table(I.PP)))[2:1],1),NA),
              c(OR,NA,NA),
              c(IC[1],NA,NA),
              c(IC[2],NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Inf. Anteparto (N)", "Inf. Anteparto (%)","Inf. Periparto (N)", "Inf. Periparto (%)","OR","Extr. Inf. IC","Extr. Sup. IC", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
############
TablaAntePeriOR.m=function(I,L,r=6){
  I.AP=I[Casos$PreP==1]
  I.PP=I[Casos$PreP==0]
  n_IAP=length(I.AP)
  n_IPP=length(I.PP)
  NA_IAP=length(I.AP[is.na(I.AP)])
  NA_IPP=length(I.PP[is.na(I.PP)])
  EE=rbind(as.vector(table(I.AP)),
           as.vector(table(I.PP)))
  pp=rep(0,dim(EE)[2])
  OR=rep(0,dim(EE)[2])
  CI1=rep(0,dim(EE)[2])
  CI2=rep(0,dim(EE)[2])
  
  for (i in 1:dim(EE)[2]){
    EE.temp=rbind(c(EE[1,i],sum(EE[1,-i])),
                  c(EE[2,i],sum(EE[2,-i])))
    FT=fisher.test(EE.temp)
    pp[i]=FT$p.value
    OR[i]=round(FT$estimate,2)
    CI1[i]=round(FT$conf.int,2)[1]
    CI2[i]=round(FT$conf.int,2)[2]
  }
  pp=p.adjust(pp,method="bonferroni")
  
  EEExt=rbind(c(as.vector(table(I.AP)),NA_IAP),
              c(round(100*as.vector(prop.table(table(I.AP))),1),NA),
              c(as.vector(table(I.PP)),NA_IPP),
              c(round(100*as.vector(prop.table(table(I.PP))),1),NA),
              c(OR,NA),
              c(CI1,NA),
              c(CI2,NA),
              c(round(pp,r),NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Inf. Anteparto (N)", "Inf. Anteparto (%)","Inf. Periparto (N)", "Inf. Periparto (%)","OR","Extr. Inf. IC","Extr. Sup. IC", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
#############


TablaAntePeri=function(I,L2,L1,r=6){
  I.AP=I[Casos$PreP==1]
  I.PP=I[Casos$PreP==0]
  n_IAP=length(I.AP)
  n_IPP=length(I.PP)
  NA_IAP=length(I.AP[is.na(I.AP)])
  NA_IPP=length(I.PP[is.na(I.PP)])
  EE=rbind(as.vector(table(I.AP))[2:1],
           as.vector(table(I.PP))[2:1])
  FT=fisher.test(EE)
  PT=prop.test(EE[,1],rowSums(EE))
  
  RA=round(PT$estimate[1]-PT$estimate[2],4)
  ICRA=round(PT$conf.int,4)
  RR=round(PT$estimate[1]/PT$estimate[2],2)
  ICRR=round(RelRisk(EE,conf.level=0.95,method="score")[2:3],2)
  
  EEExt=rbind(c(as.vector(table(I.AP))[2:1],NA_IAP),
              c(round(100*as.vector(prop.table(table(I.AP)))[2:1],1),NA),
              c(as.vector(table(I.PP))[2:1],NA_IPP),
              c(round(100*as.vector(prop.table(table(I.PP)))[2:1],1),NA),
              c(RA,NA,NA),
              c(ICRA[1],NA,NA),
              c(ICRA[2],NA,NA),
              c(RR,NA,NA),
              c(ICRR[1],NA,NA),
              c(ICRR[2],NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Inf. Anteparto (N)", "Inf. Anteparto (%)","Inf. Periparto (N)", "Inf. Periparto (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
####
TablaAntePeri.m=function(I,L,r=6){
  I.AP=I[Casos$PreP==1]
  I.PP=I[Casos$PreP==0]
  n_IAP=length(I.AP)
  n_IPP=length(I.PP)
  NA_IAP=length(I.AP[is.na(I.AP)])
  NA_IPP=length(I.PP[is.na(I.PP)])
  EE=rbind(as.vector(table(I.AP)),
           as.vector(table(I.PP)))
  pp=rep(0,dim(EE)[2])
  RA=rep(0,dim(EE)[2])
  RR=rep(0,dim(EE)[2])
  CIRA1=rep(0,dim(EE)[2])
  CIRA2=rep(0,dim(EE)[2])
  CIRR1=rep(0,dim(EE)[2])
  CIRR2=rep(0,dim(EE)[2])
  
  for (i in 1:dim(EE)[2]){
    EE.temp=rbind(c(EE[1,i],sum(EE[1,-i])),
                  c(EE[2,i],sum(EE[2,-i])))
    FT=fisher.test(EE.temp)
    PT=prop.test(EE.temp[,1],rowSums(EE.temp))
    pp[i]=FT$p.value
    RA[i]=round(PT$estimate[1]-PT$estimate[2],4)
    RR[i]=round(PT$estimate[1]/PT$estimate[2],3)
    CIRA1[i]=round(PT$conf.int,4)[1]
    CIRA2[i]=round(PT$conf.int,4)[2]
    CIRR1[i]=round(RelRisk(EE.temp,conf.level=0.95,method="wald"),3)[2]
    CIRR2[i]=round(RelRisk(EE.temp,conf.level=0.95,method="wald"),3)[3]
  }
  pp=p.adjust(pp,method="bonferroni")
  
  EEExt=rbind(c(as.vector(table(I.AP)),NA_IAP),
              c(round(100*as.vector(prop.table(table(I.AP))),1),NA),
              c(as.vector(table(I.PP)),NA_IPP),
              c(round(100*as.vector(prop.table(table(I.PP))),1),NA),
              c(RA,NA),
              c(CIRA1,NA),
              c(CIRA2,NA),
              c(RR,NA),
              c(CIRR1,NA),
              c(CIRR2,NA),
              c(round(pp,r),NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Inf. Anteparto (N)", "Inf. Anteparto (%)","Inf. Periparto (N)", "Inf. Periparto (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
#############
###########
TablaAntePeri.r=function(I.AP,I.PP,L2,L1,r=6){
  n_IAP=length(I.AP)
  n_IPP=length(I.PP)
  NA_IAP=length(I.AP[is.na(I.AP)])
  NA_IPP=length(I.PP[is.na(I.PP)])
  EE=rbind(as.vector(table(I.AP))[2:1],
           as.vector(table(I.PP))[2:1])
  FT=fisher.test(EE)
  PT=prop.test(EE[,1],rowSums(EE))
  
  RA=round(PT$estimate[1]-PT$estimate[2],4)
  ICRA=round(PT$conf.int,4)
  RR=round(PT$estimate[1]/PT$estimate[2],2)
  ICRR=round(RelRisk(EE,conf.level=0.95,method="score")[2:3],2)
  
  EEExt=rbind(c(as.vector(table(I.AP))[2:1],NA_IAP),
              c(round(100*as.vector(prop.table(table(I.AP)))[2:1],1),NA),
              c(as.vector(table(I.PP))[2:1],NA_IPP),
              c(round(100*as.vector(prop.table(table(I.PP)))[2:1],1),NA),
              c(RA,NA,NA),
              c(ICRA[1],NA,NA),
              c(ICRA[2],NA,NA),
              c(RR,NA,NA),
              c(ICRR[1],NA,NA),
              c(ICRR[2],NA,NA),
              c(round(FT$p.value,r),NA,NA)
  )
  
  colnames(EEExt)=c(L2,L1,"Datos perdidos")
  rownames(EEExt)=c("Inf. Anteparto (N)", "Inf. Anteparto (%)","Inf. Periparto (N)", "Inf. Periparto (%)","RA","Extr. Inf. IC RA","Extr. Sup. IC RA","RR","Extr. Inf. IC RA","Extr. Sup. IC RR", "p-valor")
  
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
#####
#####

TablaAntePeriCasos=function(CasosX,I,L,SintX=Sint.Tot){
  I.AP.A=I[CasosX$PreP==1& SintX=="Asintomática"]
  I.AP.L=I[CasosX$PreP==1& SintX=="Leve"]
  I.AP.G=I[CasosX$PreP==1& SintX=="Grave"]
  I.PP.A=I[CasosX$PreP==0& SintX=="Asintomática"]
  I.PP.L=I[CasosX$PreP==0& SintX=="Leve"]
  I.PP.G=I[CasosX$PreP==0& SintX=="Asintomática"]
  
  NA_I.AP.A=length(I.AP.A[is.na(I.AP.A)])
  NA_I.AP.L=length(I.AP.L[is.na(I.AP.L)])
  NA_I.AP.G=length(I.AP.G[is.na(I.AP.G)])
  NA_I.PP.A=length(I.PP.A[is.na(I.PP.A)])
  NA_I.PP.L=length(I.PP.L[is.na(I.PP.L)])
  NA_I.PP.G=length(I.PP.G[is.na(I.PP.G)])
  
  
  
  EEExt=rbind(c(as.vector(table(I.AP.A))[2:1],NA_I.AP.A),
              c(round(100*as.vector(prop.table(table(I.AP.A)))[2:1],1),NA),
              c(as.vector(table(I.AP.L))[2:1],NA_I.AP.L),
              c(round(100*as.vector(prop.table(table(I.AP.L)))[2:1],1),NA),
              c(as.vector(table(I.AP.G))[2:1],NA_I.AP.G),
              c(round(100*as.vector(prop.table(table(I.AP.G)))[2:1],1),NA),
              c(as.vector(table(I.PP.A))[2:1],NA_I.PP.A),
              c(round(100*as.vector(prop.table(table(I.PP.A)))[2:1],1),NA),
              c(as.vector(table(I.PP.L))[2:1],NA_I.PP.L),
              c(round(100*as.vector(prop.table(table(I.PP.L)))[2:1],1),NA),
              c(as.vector(table(I.PP.G))[2:1],NA_I.PP.G),
              c(round(100*as.vector(prop.table(table(I.PP.G)))[2:1],1),NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Inf. AP. Asintomáticas (N)", "Inf. AP. Asintomáticas (%)",
                    "Inf. AP. Leves (N)", "Inf. AP. Leves (%)",
                    "Inf. AP. Graves (N)", "Inf. AP. Graves (%)",
                    "Inf. PP. Asintomáticas (N)", "Inf. PP. Asintomáticas (%)",
                    "Inf. PP. Leves (N)", "Inf. PP. Leves (%)",
                    "Inf. PP. Graves (N)", "Inf. PP. Graves (%)")
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}
#########
########
TablaAntePeriCasos.m=function(CasosX,I,L,SintX=Sint.Tot){
  I.AP.A=I[CasosX$PreP==1& Sint=="Asintomática"]
  I.AP.L=I[CasosX$PreP==1& SintX=="Leve"]
  I.AP.G=I[CasosX$PreP==1& SintX=="Grave"]
  I.PP.A=I[CasosX$PreP==0& SintX=="Asintomática"]
  I.PP.L=I[CasosX$PreP==0& SintX=="Leve"]
  I.PP.G=I[CasosX$PreP==0& SintX=="Asintomática"]
  
  NA_I.AP.A=length(I.AP.A[is.na(I.AP.A)])
  NA_I.AP.L=length(I.AP.L[is.na(I.AP.L)])
  NA_I.AP.G=length(I.AP.G[is.na(I.AP.G)])
  NA_I.PP.A=length(I.PP.A[is.na(I.PP.A)])
  NA_I.PP.L=length(I.PP.L[is.na(I.PP.L)])
  NA_I.PP.G=length(I.PP.G[is.na(I.PP.G)])
  
  
  
  EEExt=rbind(c(as.vector(table(I.AP.A)),NA_I.AP.A),
              c(round(100*as.vector(prop.table(table(I.AP.A))),1),NA),
              c(as.vector(table(I.AP.L)),NA_I.AP.L),
              c(round(100*as.vector(prop.table(table(I.AP.L))),1),NA),
              c(as.vector(table(I.AP.G)),NA_I.AP.G),
              c(round(100*as.vector(prop.table(table(I.AP.G))),1),NA),
              c(as.vector(table(I.PP.A)),NA_I.PP.A),
              c(round(100*as.vector(prop.table(table(I.PP.A))),1),NA),
              c(as.vector(table(I.PP.L)),NA_I.PP.L),
              c(round(100*as.vector(prop.table(table(I.PP.L))),1),NA),
              c(as.vector(table(I.PP.G)),NA_I.PP.G),
              c(round(100*as.vector(prop.table(table(I.PP.G))),1),NA)
  )
  
  colnames(EEExt)=c(L,"Datos perdidos")
  rownames(EEExt)=c("Inf. AP. Asintomáticas (N)", "Inf. AP. Asintomáticas (%)",
                    "Inf. AP. Leves (N)", "Inf. AP. Leves (%)",
                    "Inf. AP. Graves (N)", "Inf. AP. Graves (%)",
                    "Inf. PP. Asintomáticas (N)", "Inf. PP. Asintomáticas (%)",
                    "Inf. PP. Leves (N)", "Inf. PP. Leves (%)",
                    "Inf. PP. Graves (N)", "Inf. PP. Graves (%)")
  
  t(EEExt) %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}

########
# Tablas 2d Robson

TablaRobsonGlobal=function(TR1,TR2,x,y=x){
  EE1=t(TR1)[2:1,2:11]
  EE2=t(TR2)[2:1,2:11]
  OR=rep(0,10)
  pp=rep(0,10)
  for (i in 1:10){
    FT=fisher.test(cbind(EE1[,i],EE2[,i]))
    OR[i]=FT$estimate
    pp[i]=FT$p.value       }
  
  dt=data.frame(1:10,colSums(EE1),EE1[1,],round(100*prop.table(EE1,margin=2)[1,],2),colSums(EE2),EE2[1,],round(100*prop.table(EE2,margin=2)[1,],2),round(OR,2),round(pp,3))
  names(dt)=c("Grupo Robson", paste("Casos",x), "Cesáreas", "Porcentaje",  paste("Controles",y), "Cesáreas", "Porcentaje","OR","p-valor test de Fisher")
  dt %>%
    kbl() %>%
    kable_styling()%>%     
    scroll_box(width="100%", box_css="border: 0px;")}



# Barplots descripción global sintomatología

Barplot.DGS=function(I,L1,L2){
  taula=table(data.frame(Factor=I,Síntomas=Sint))[c(2,1) ,1:3]
  Síntomas=ordered(rep(c("Asintomática", "Leve", "Grave"), each=2),levels=c("Asintomática", "Leve", "Grave"))
  Grupo=ordered(rep(c(L1, L2) , 3),levels=c(L1, L2))
  valor=as.vector(prop.table(taula, margin=1))
  data <- data.frame(Grupo,Síntomas,valor)
  
  ggplot(data, aes(fill=Síntomas, y=valor, x=Grupo)) + 
    geom_bar(position="dodge", stat="identity")+
    ylab("")+
    xlab("")+  
    scale_fill_brewer(palette = "Dark2") 
}

Barplot.DGSalt=function(I,S,L1,L2){
  taula=table(data.frame(Factor=I,Síntomas=S))[c(2,1) ,1:3]
  Síntomas=ordered(rep(c("Asintomática", "Leve", "Grave"), each=2),levels=c("Asintomática", "Leve", "Grave"))
  Grupo=ordered(rep(c(L1, L2) , 3),levels=c(L1, L2))
  valor=as.vector(prop.table(taula, margin=1))
  data <- data.frame(Grupo,Síntomas,valor)
  
  ggplot(data, aes(fill=Síntomas, y=valor, x=Grupo)) + 
    geom_bar(position="dodge", stat="identity")+
    ylab("")+
    xlab("")+  
    scale_fill_brewer(palette = "Dark2") 
}


Barplot.DGSW=function(I,L1,L2){
  taula=table(data.frame(Factor=I,Síntomas=Sint))[c(2,1),]
  Síntomas=ordered(rep(c("Asintomática", "Leve", "Grave"), each=2),levels=c("Asintomática", "Leve", "Grave"))
  Grupo=ordered(rep(c(L1,L2) , 3),levels=c(L1,L2))
  valor=as.vector(prop.table(taula, margin=1))
  data <- data.frame(Grupo,Síntomas,valor)
  
  ggplot(data, aes(fill=Síntomas, y=valor, x=Grupo)) + 
    geom_bar(position="dodge", stat="identity")+
    ylab("")+
    xlab("")+  
    scale_fill_brewer(palette = "Dark2") 
}




ComparOlasCC=function(I1,I2,NI1,NI2,L1,L2,r=6){
  Tot1w=c(I1,NI1)
  Tot2w=c(I2,NI2)
  NA_I1=length(I1[is.na(I1)])
  NA_NI1=length(NI1[is.na(NI1)])
  NA_I2=length(I2[is.na(I2)])
  NA_NI2=length(NI2[is.na(NI2)])
  EE=cbind(as.vector(table(I1))[2:1],
           as.vector(table(NI1))[2:1],
           as.vector(table(I2))[2:1],
           as.vector(table(NI2))[2:1])
  
  OR1=fisher.test(EE[,1:2])$estimate
  OR2=fisher.test(EE[,3:4])$estimate
  QOR=OR1/OR2
  pp=round(OR.OR(EE)[1],r)
  IC=OR.OR(EE)[2:3]
  
  EETot=cbind(as.vector(table(Tot1w)),as.vector(table(Tot2w)))
  EEI=EE[,c(1,3)]
  EENI=EE[,c(2,4)]
  FTGlob=c(fisher.test(EETot)$p.value,fisher.test(EEI)$p.value,fisher.test(EENI)$p.value)
  
  EEExt=cbind(c(as.vector(table(I1))[2:1],NA_I1),
              c(round(100*as.vector(prop.table(table(I1)))[2:1],1),NA),
              c(as.vector(table(NI1))[2:1],NA_NI1),
              c(round(100*as.vector(prop.table(table(NI1)))[2:1],1),NA),
              c(as.vector(table(I2))[2:1],NA_I2),
              c(round(100*as.vector(prop.table(table(I2)))[2:1],1),NA),
              c(as.vector(table(NI2))[2:1],NA_NI2),
              c(round(100*as.vector(prop.table(table(NI2)))[2:1],1),NA),
              c(round(QOR,4),NA,NA),
              c(round(IC[1],4),NA,NA),
              c(round(IC[2],4),NA,NA),
              c(pp,NA,NA),
              c(round(FTGlob[1],r),NA,NA),
              c(round(FTGlob[2],r),NA,NA),
              c(round(FTGlob[3],r),NA,NA)
  )
  colnames(EEExt)=c(ColumnesExt,"Cociente de OR", "Extr. inf. IC",  "Extr. sup. IC", "p-valor igualdad OR","p-valor global olas", "p-valor casos","p-valor controles")
  
  rownames(EEExt)=c(L1, L2,"Datos perdidos")
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}



ComparOlasCCm=function(I1,I2,NI1,NI2,L,r=6){
  Tot1w=c(I1,NI1)
  Tot2w=c(I2,NI2)
  NA_I1=length(I1[is.na(I1)])
  NA_NI1=length(NI1[is.na(NI1)])
  NA_I2=length(I2[is.na(I2)])
  NA_NI2=length(NI2[is.na(NI2)])
  
  EE=cbind(as.vector(table(I1)),
           as.vector(table(NI1)),
           as.vector(table(I2)),
           as.vector(table(NI2)))
  EETot=cbind(as.vector(table(Tot1w)),as.vector(table(Tot2w)))
  
  
  OR1=rep(NA,dim(EE)[1])
  OR2=rep(NA,dim(EE)[1])
  QOR=rep(NA,dim(EE)[1])
  pp=rep(NA,dim(EE)[1])
  IC1=rep(NA,dim(EE)[1])
  IC2=rep(NA,dim(EE)[1])
  FTGlob1=rep(NA,dim(EE)[1])
  FTGlob2=rep(NA,dim(EE)[1])
  FTGlob3=rep(NA,dim(EE)[1])
  
  for (i in 1: dim(EE)[1]){
    EEt=rbind(EE[i,],colSums(EE[-i,]))
    EETott=rbind(EETot[i,],colSums(EETot[-i,])) 
    QOR[i]=fisher.test(EEt[,1:2])$estimate/fisher.test(EEt[,3:4])$estimate
    pp[i]=OR.OR(EEt)[1]
    IC1[i]=OR.OR(EEt)[2]
    IC2[i]=OR.OR(EEt)[3]
    FTGlob1[i]=fisher.test(EETott)$p.value
    FTGlob2[i]=fisher.test(EEt[,c(1,3)])$p.value
    FTGlob3[i]=fisher.test(EEt[,c(2,4)])$p.value
  }
  pp=round(p.adjust(pp,method="bonferroni"),r)
  FTGlob1=round(p.adjust(FTGlob1,method="bonferroni"),r)
  FTGlob2=round(p.adjust(FTGlob2,method="bonferroni"),r)
  FTGlob3=round(p.adjust(FTGlob3,method="bonferroni"),r)
  
  
  EEExt=cbind(c(as.vector(table(I1)),NA_I1),
              c(round(100*as.vector(prop.table(table(I1))),1),NA),
              c(as.vector(table(NI1)),NA_NI1),
              c(round(100*as.vector(prop.table(table(NI1))),1),NA),
              c(as.vector(table(I2)),NA_I2),
              c(round(100*as.vector(prop.table(table(I2))),1),NA),
              c(as.vector(table(NI2)),NA_NI2),
              c(round(100*as.vector(prop.table(table(NI2))),1),NA),
              c(round(QOR,2),NA),
              c(round(IC1,2),NA),
              c(round(IC2,2),NA),
              c(pp,NA),
              c(FTGlob1,NA),
              c(FTGlob2,NA),
              c(FTGlob3,NA)
  )
  colnames(EEExt)=c(ColumnesExt,"Cociente de OR", "Extr. inf. IC",  "Extr. sup. IC", "p-valor igualdad OR","p-valor global olas", "p-valor casos","p-valor controles")
  
  rownames(EEExt)=c(L,"Datos perdidos")
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}





ComparOlasRR=function(I1,I2,NI1,NI2,L1,L2,r=6){
  Tot1w=c(I1,NI1)
  Tot2w=c(I2,NI2)
  NA_I1=length(I1[is.na(I1)])
  NA_NI1=length(NI1[is.na(NI1)])
  NA_I2=length(I2[is.na(I2)])
  NA_NI2=length(NI2[is.na(NI2)])
  EE=cbind(as.vector(table(I1))[2:1],
           as.vector(table(NI1))[2:1],
           as.vector(table(I2))[2:1],
           as.vector(table(NI2))[2:1])
  
  e=EE[1,]
  n=colSums(EE)
  p=e/n
  
  RA=round((p[1]-p[2])-(p[3]-p[4]),4)
  RR=round((p[1]/p[2])/(p[3]/p[4]),4)
  p.RARA=RA.RA(EE)[1]
  IC.RARA=RA.RA(EE)[2:3]
  p.RRRR=RR.RR(EE)[2]
  IC.RRRR=RR.RR(EE)[3:4]
  
  
  EETot=cbind(as.vector(table(Tot1w)),as.vector(table(Tot2w)))
  EEI=EE[,c(1,3)]
  EENI=EE[,c(2,4)]
  FTGlob=c(fisher.test(EETot)$p.value,fisher.test(EEI)$p.value,fisher.test(EENI)$p.value)
  
  EEExt=cbind(c(as.vector(table(I1))[2:1],NA_I1),
              c(round(100*as.vector(prop.table(table(I1)))[2:1],1),NA),
              c(as.vector(table(NI1))[2:1],NA_NI1),
              c(round(100*as.vector(prop.table(table(NI1)))[2:1],1),NA),
              c(as.vector(table(I2))[2:1],NA_I2),
              c(round(100*as.vector(prop.table(table(I2)))[2:1],1),NA),
              c(as.vector(table(NI2))[2:1],NA_NI2),
              c(round(100*as.vector(prop.table(table(NI2)))[2:1],1),NA),
              c(round(RA,4),NA,NA),
              c(round(IC.RARA[1],4),NA,NA),
              c(round(IC.RARA[2],4),NA,NA),
              c(round(p.RARA,r),NA,NA),
              c(round(RR,4),NA,NA),
              c(round(IC.RRRR[1],4),NA,NA),
              c(round(IC.RRRR[2],4),NA,NA),
              c(round(p.RRRR,r),NA,NA),
              c(round(FTGlob[1],r),NA,NA),
              c(round(FTGlob[2],r),NA,NA),
              c(round(FTGlob[3],r),NA,NA)
  )
  colnames(EEExt)=c(ColumnesExt,"Diferencia RA", "Extr. inf. IC RA",  "Extr. sup. IC RA", "p-valor igualdad RA","Cociente RR", "Extr. inf. IC RR",  "Extr. sup. IC RR", "p-valor igualdad RR","p-valor global olas", "p-valor casos","p-valor controles")
  
  rownames(EEExt)=c(L1, L2,"Datos perdidos")
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}




ComparOlasRRwe=function(I1,I2,NI1,NI2,L1,L2,r=6){
  Tot1w=c(I1,NI1)
  Tot2w=c(I2,NI2)
  NA_I1=length(I1[is.na(I1)])
  NA_NI1=length(NI1[is.na(NI1)])
  NA_I2=length(I2[is.na(I2)])
  NA_NI2=length(NI2[is.na(NI2)])
  EE=cbind(as.vector(table(I1))[2:1],
           as.vector(table(NI1))[2:1],
           as.vector(table(I2))[2:1],
           as.vector(table(NI2))[2:1])
  
  e=EE[1,]
  n=colSums(EE)
  p=e/n
  
  RA=round((p[1]-p[2])-(p[3]-p[4]),4)
  RR=round((p[1]/p[2])/(p[3]/p[4]),4)
  p.RARA=RA.RA(EE)[1]
  IC.RARA=RA.RA(EE)[2:3]
  p.RRRR=RR.RR(EE)[2]
  IC.RRRR=RR.RR(EE)[3:4]
  
  
  EETot=cbind(as.vector(table(Tot1w)),as.vector(table(Tot2w)))
  EEI=EE[,c(1,3)]
  EENI=EE[,c(2,4)]
  FTGlob=c(fisher.test(EETot,workspace=2e8)$p.value,fisher.test(EEI,workspace=2e8)$p.value,fisher.test(EENI,workspace=2e8)$p.value)
  
  EEExt=cbind(c(as.vector(table(I1))[2:1],NA_I1),
              c(round(100*as.vector(prop.table(table(I1)))[2:1],1),NA),
              c(as.vector(table(NI1))[2:1],NA_NI1),
              c(round(100*as.vector(prop.table(table(NI1)))[2:1],1),NA),
              c(as.vector(table(I2))[2:1],NA_I2),
              c(round(100*as.vector(prop.table(table(I2)))[2:1],1),NA),
              c(as.vector(table(NI2))[2:1],NA_NI2),
              c(round(100*as.vector(prop.table(table(NI2)))[2:1],1),NA),
              c(round(RA,4),NA,NA),
              c(round(IC.RARA[1],4),NA,NA),
              c(round(IC.RARA[2],4),NA,NA),
              c(round(p.RARA,r),NA,NA),
              c(round(RR,4),NA,NA),
              c(round(IC.RRRR[1],4),NA,NA),
              c(round(IC.RRRR[2],4),NA,NA),
              c(round(p.RRRR,r),NA,NA),
              c(round(FTGlob[1],r),NA,NA),
              c(round(FTGlob[2],r),NA,NA),
              c(round(FTGlob[3],r),NA,NA)
  )
  colnames(EEExt)=c(ColumnesExt,"Diferencia RA", "Extr. inf. IC RA",  "Extr. sup. IC RA", "p-valor igualdad RA","Cociente RR", "Extr. inf. IC RR",  "Extr. sup. IC RR", "p-valor igualdad RR","p-valor global olas", "p-valor casos","p-valor controles")
  
  rownames(EEExt)=c(L1, L2,"Datos perdidos")
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}



ComparOlasRR.m=function(I1,I2,NI1,NI2,L,r=6){
  Tot1w=c(I1,NI1)
  Tot2w=c(I2,NI2)
  NA_I1=length(I1[is.na(I1)])
  NA_NI1=length(NI1[is.na(NI1)])
  NA_I2=length(I2[is.na(I2)])
  NA_NI2=length(NI2[is.na(NI2)])
  
  
  EE=cbind(as.vector(table(I1)),
           as.vector(table(NI1)),
           as.vector(table(I2)),
           as.vector(table(NI2)))
  
  EETot=cbind(as.vector(table(Tot1w)),as.vector(table(Tot2w)))
  
  RA=rep(NA,dim(EE)[1])
  RR=rep(NA,dim(EE)[1])
  p.RARA=rep(NA,dim(EE)[1])
  IC.RARA.1=rep(NA,dim(EE)[1])
  IC.RARA.2=rep(NA,dim(EE)[1])
  p.RRRR=rep(NA,dim(EE)[1])
  IC.RRRR.1=rep(NA,dim(EE)[1])
  IC.RRRR.2=rep(NA,dim(EE)[1])
  
  FTGlob1=rep(NA,dim(EE)[1])
  FTGlob2=rep(NA,dim(EE)[1])
  FTGlob3=rep(NA,dim(EE)[1])
  
  for (i in 1: dim(EE)[1]){
    EEt=rbind(EE[i,],colSums(EE[-i,]))
    EETott=rbind(EETot[i,],colSums(EETot[-i,])) 
    e=EEt[1,]
    n=colSums(EEt)
    p=e/n
    
    RA[i]=round((p[1]-p[2])-(p[3]-p[4]),4)
    RR[i]=round((p[1]/p[2])/(p[3]/p[4]),4)
    p.RARA[i]=RA.RA(EEt)[1]
    IC.RARA.1[i]=RA.RA(EEt)[2]
    IC.RARA.2[i]=RA.RA(EEt)[3]
    p.RRRR[i]=RR.RR(EEt)[2]
    IC.RRRR.1[i]=RR.RR(EEt)[3]
    IC.RRRR.2[i]=RR.RR(EEt)[4]
    
    FTGlob1[i]=fisher.test(EETott)$p.value
    FTGlob2[i]=fisher.test(EEt[,c(1,3)])$p.value
    FTGlob3[i]=fisher.test(EEt[,c(2,4)])$p.value
  }
  p.RARA=round(p.adjust(p.RARA,method="bonferroni"),r)
  p.RRRR=round(p.adjust(p.RRRR,method="bonferroni"),r)
  FTGlob1=round(p.adjust(FTGlob1,method="bonferroni"),r)
  FTGlob2=round(p.adjust(FTGlob2,method="bonferroni"),r)
  FTGlob3=round(p.adjust(FTGlob3,method="bonferroni"),r)
  
  
  
  
  EEExt=cbind(c(as.vector(table(I1)),NA_I1),
              c(round(100*as.vector(prop.table(table(I1))),1),NA),
              c(as.vector(table(NI1)),NA_NI1),
              c(round(100*as.vector(prop.table(table(NI1))),1),NA),
              c(as.vector(table(I2)),NA_I2),
              c(round(100*as.vector(prop.table(table(I2))),1),NA),
              c(as.vector(table(NI2)),NA_NI2),
              c(round(100*as.vector(prop.table(table(NI2))),1),NA),
              c(round(RA,2),NA),
              c(round(IC.RARA.1,2),NA),
              c(round(IC.RARA.2,2),NA),
              c(p.RARA,NA),
              c(round(RR,2),NA),
              c(round(IC.RRRR.1,2),NA),
              c(round(IC.RRRR.2,2),NA),
              c(p.RRRR,NA),
              c(FTGlob1,NA),
              c(FTGlob2,NA),
              c(FTGlob3,NA)
  )
  colnames(EEExt)=c(ColumnesExt,"Diferencia RA", "Extr. inf. IC RA", 
                    "Extr. sup. IC RA", "p-valor igualdad RA","Cociente RR", "Extr. inf. IC RR",  "Extr. sup. IC RR", "p-valor igualdad RR","p-valor global olas", "p-valor casos","p-valor controles")
  rownames(EEExt)=c(L,"Datos perdidos")
  
  EEExt %>%
    kbl() %>%
    kable_styling() %>%     
    scroll_box(width="100%", box_css="border: 0px;")  
}



GrComparOlasCC=function(I1,I2,NI1,NI2,L1){
  base=c("1a ola","2a ola")
  df=data.frame(
    Factor=c(rep("Infectadas 1a ola",length(I1)), rep("No infectadas 1a ola",length(NI1)),rep("Infectadas 2a ola",length(I2)), rep("No infectadas 2a ola",length(NI2))),
    Grupo=c(I1,NI1,I2,NI2 )
  )
  Grupo=ordered(rep(c("Infectadas","No infectadas"), 2),levels=c("Infectadas","No infectadas"))
  Ola=ordered(rep(base , each=2),levels=base)
  valor=as.vector(100*prop.table(table(df)[c(1,3,2,4),2:1], margin=1))[1:4]
  data=data.frame(Grupo,Ola,valor)
  data %>%
    ggplot(aes(fill=Grupo, y=valor, x=Ola)) + 
    geom_bar(position="dodge", stat="identity")+
    xlab("")+
    ylab(paste(L1,"(%)"))+
    scale_fill_hue(c=100)+labs(fill="")}


