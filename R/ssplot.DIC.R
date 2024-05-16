

ssplot.DIC=function(D,design=2,quali=c(TRUE,TRUE,TRUE),sig=0.05,verbose=TRUE){

FA=D[,1]
FB=D[,2]
FC=D[,3]
Rep=D[,4]
Y=D[,5]
vari=colnames(D)[5]



a=length(unique(FA))
b=length(unique(FB))
c=length(unique(FC))


bloco=sq(Rep,Y)
fa=sq(FA,Y)
parc1=sq(paste(FA,Rep),Y)
erroA=parc1-fa

fb=sq(FB,Y)
fab=sqI2(FA,FB,Y)
parc2=sq(paste(FA,FB,Rep),Y)
erroB=parc2-parc1-fb-fab

fc=sq(FC,Y)
fac=sqI2(FA,FC,Y)
fbc=sqI2(FB,FC,Y)
fabc=sqI3(FA,FB,FC,Y)
total=sq(1:length(Y),Y)
erroc=total-parc2-fc-fac-fbc-fabc



anova=rbind(bloco,fa,erroA,fb,fab,erroB,fc,fac,fbc,fabc,erroc)
rownames(anova)=c( "bloco" ,"A",    "Residual A", "B",    "AxB" ,  "Residual B" ,"C" ,
                   "AxC",   "BxC"  , "AxBxC",  "Residual c")



QM=anova[,2]/anova[,1]
anova=cbind(anova,QM=QM)
Fvalue=c(NA,QM[2]/QM[3],NA,QM[4]/QM[6],QM[5]/QM[6],NA,QM[7]/QM[11],QM[8]/QM[11],QM[9]/QM[11],QM[10]/QM[11],NA)
Fvalue=round(Fvalue,5)
pvalue=c(NA,
        1-pf(QM[2]/QM[3],anova[2,1],anova[3,1]),
         NA,
        1-pf(QM[4]/QM[6],anova[4,1],anova[6,1]),
        1-pf(QM[5]/QM[6],anova[5,1],anova[6,1]),
         NA,
        1-pf(QM[7]/QM[11],anova[7,1],anova[11,1]),
        1-pf(QM[8]/QM[11],anova[8,1],anova[11,1]),
        1-pf(QM[9]/QM[11],anova[9,1],anova[11,1]),
        1-pf(QM[10]/QM[11],anova[10,1],anova[11,1]),
         NA)

pvalue2=round(pvalue,5)
pvalue2[pvalue2<0.0001]="<0.0001"

anova=data.frame(anova,Fvalue,pvalue=pvalue2)
anova[is.na(anova)]=""
# GlFa
# QMFA
#
# GlFb
# QMFb
#
# GlFc
# QMFc
anova=anova[-1,]

glea=anova[2,1]
gleb=anova[5,1]
glec=anova[10,1]

qmea=anova[2,3]
qmeb=anova[5,3]
qmec=anova[10,3]

CV1=round(100*sqrt(qmea)/mean(Y),5)
CV2=round(100*sqrt(qmeb)/mean(Y),5)
CV3=round(100*sqrt(qmec)/mean(Y),5)


if(verbose==TRUE){ cat("#################################################################", "\n")
  cat(paste("Analysis of the variable (Analise da variavel) -> ",vari), "\n")
  cat(" ", "\n")
   print(anova)
   cat(" ", "\n")
   cat("Residual coefficients of variation (Coeficientes de variacao residual)", "\n")
   print(paste("CV1(%) =",CV1))
   print(paste("CV2(%) =",CV2))
   print(paste("CV3(%) =",CV3))
  cat("#################################################################", "\n")}





if(verbose==TRUE){ cat("#################################################################", "\n")
  cat("Study of the main effects (Estudo dos efeitos principais)", "\n")

  cat("#################################################################", "\n")}



#Teste Fator A
if(verbose==TRUE){ cat(".................................................................", "\n")}
if(quali[1]==TRUE){
  AnalA=ComparacaoMedias(Y,FA,DFerror = glea,MSerror = qmea,alpha = sig)
}






if(verbose==TRUE){
  cat("Comparison of levels of factor A (Comparacao dos niveis do fator A)", "\n")
  if(quali[1]==FALSE){

    fa=sq(FA,Y)
    AnalA=RegressaoPolinomial(resp=Y, trat=FA, glres=glea, SQres=qmea*glea, gltrat=fa[1], SQtrat=fa[2],verbose=T)
    AnalA=""
    }
  print(AnalA)
  }


#Teste Fator B
if(verbose==TRUE){ cat(".................................................................", "\n")}
if(quali[2]==TRUE){
  AnalB=ComparacaoMedias(Y,FB,DFerror = gleb,MSerror = qmeb,alpha = sig)
}



if(verbose==TRUE){
  cat("Comparison of levels of factor B (Comparacao dos niveis do fator B)", "\n")
  if(quali[2]==FALSE){
    fa=sq(FB,Y)
    AnalB=RegressaoPolinomial(resp=Y, trat=FB, glres=gleb, SQres=qmeb*gleb, gltrat=fa[1], SQtrat=fa[2],verbose=T)
    AnalB=""
     }
  print(AnalB)
}

#Teste Fator C
if(verbose==TRUE){ cat(".................................................................", "\n")}
if(quali[3]==TRUE){
  AnalC=ComparacaoMedias(Y,FC,DFerror = glec,MSerror = qmec,alpha = sig)
}



if(verbose==TRUE){
  cat("Comparison of levels of factor C (Comparacao dos niveis do fator C)", "\n")
  if(quali[3]==FALSE){
    fa=sq(FC,Y)
    AnalC=RegressaoPolinomial(resp=Y, trat=FC, glres=glec, SQres=qmec*glec, gltrat=fa[1], SQtrat=fa[2],verbose=T)
    AnalC=""
    }
  print(AnalC)
}



if(verbose==TRUE){ cat("#################################################################", "\n")
  cat("Study of double interactions (Estudo das interacoes duplas)", "\n")
  cat("#################################################################", "\n")}


#Desdobramento A/B
glem=Satterthwaite(QMA = qmea,QMB = qmeb,DFa = glea,DFb = gleb,b = b)
qmem=(qmea+(b-1)*qmeb)/b

for(i in unique(FB)){
  id=FB==i


  if(quali[1]==TRUE){
  AnalAB.=ComparacaoMedias(Y[id],FA[id],DFerror = glem,MSerror = qmem,alpha = sig)
  }


  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator A dentro do nivel",i,"do fator B"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor A within level",i,"of factor B"), "\n")}

    if(verbose==TRUE){
    if(quali[1]==FALSE){
      f=sq(FA[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FA[id], glres=glem, SQres=qmem*glem, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""
       }

    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}


#Desdobramento B/A
for(i in unique(FA)){
  id=FA==i


  if(quali[2]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FB[id],DFerror = gleb,MSerror = qmeb,alpha = sig)
  }








  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator B dentro do nivel",i,"do fator A"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor B within level",i,"of factor A"), "\n")}

    if(verbose==TRUE){
    if(quali[2]==FALSE){
      f=sq(FB[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FB[id], glres=gleb, SQres=qmeb*gleb, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""}
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}


#Desdobramento A/C
glem=Satterthwaite(QMA = qmea,QMB = qmec,DFa = glea,DFb = glec,b = c)
qmem=(qmea+(c-1)*qmec)/c
for(i in unique(FC)){
  id=FC==i




  if(quali[1]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FA[id],DFerror = glem,MSerror = qmem,alpha = sig)
  }






  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator A dentro do nivel",i,"do fator C"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor A within level",i,"of factor C"), "\n")}

    if(verbose==TRUE){
    if(quali[1]==FALSE){
      f=sq(FA[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FA[id], glres=glem, SQres=qmem*glem, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""
       }
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}

#Desdobramento C/A
for(i in unique(FA)){
  id=FA==i




  if(quali[3]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FC[id],DFerror = glec,MSerror = qmec,alpha = sig)
  }



  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator C dentro do nivel",i,"do fator A"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor C within level",i,"of factor A"), "\n")}

    if(verbose==TRUE){
    if(quali[3]==FALSE){
      f=sq(FC[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FC[id], glres=glec, SQres=qmec*glec, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""
      }
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}




#Desdobramento B/C
glem=Satterthwaite(QMA = qmeb,QMB = qmec,DFa = gleb,DFb = glec,b = c)
qmem=(qmeb+(c-1)*qmec)/c
for(i in unique(FC)){
  id=FC==i



  if(quali[2]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FB[id],DFerror = glem,MSerror = qmem,alpha = sig)
  }


  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator B dentro do nivel",i,"do fator C"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor B within level",i,"of factor C"), "\n")}

    if(verbose==TRUE){
    if(quali[2]==FALSE){
      f=sq(FB[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FB[id], glres=glem, SQres=qmem*glem, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.="" }
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}

#Desdobramento C/B
for(i in unique(FB)){
  id=FB==i




  if(quali[3]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FC[id],DFerror = glec,MSerror = qmec,alpha = sig)
  }


  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator C dentro do nivel",i,"do fator B"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor C within level",i,"of factor B"), "\n")}

    if(verbose==TRUE){
    if(quali[3]==FALSE){
      f=sq(FC[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FC[id], glres=glec, SQres=qmec*glec, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""}
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}



########################################
### Interacao tripla

if(verbose==TRUE){ cat("#################################################################", "\n")
  cat("Study of triple interactions (Estudo das interacoes triplas)", "\n")
  cat("#################################################################", "\n")}

#Desdobramento A/BC
glem=Satterthwaite(QMA = qmeb,QMB = qmec,DFa = gleb,DFb = glec,b = c)
qmem=(qmeb+(c-1)*qmec)/c

glem2=Satterthwaite(QMA = qmea,QMB = glem,DFa = glea,DFb = glem,b =b*c )
qmem2=(qmea+(b*c-1)*qmem)/(b*c)


for(i in unique(FB)){
  for(j in unique(FC)){
  id=(FB==i)&(FC==j)



  if(quali[1]==TRUE){
    AnalAB.=ComparacaoMedias(Y[id],FA[id],DFerror = glem2,MSerror = qmem2,alpha = sig)
    }


  if(verbose==TRUE){ cat(".................................................................", "\n")}
  if(verbose==TRUE){ cat(paste("Desdobramento do fator A dentro do nivel",i,"do fator B e",j,"do fator C"), "\n")}
  if(verbose==TRUE){ cat(paste("Splitting of factor A within level",i,"of factor B and",j, "of factor C"), "\n")}

    if(verbose==TRUE){
    if(quali[1]==FALSE){
      f=sq(FA[id],Y[id])
      AnalAB.=RegressaoPolinomial(Y[id],FA[id], glres=glem2, SQres=qmem2*glem2, gltrat=f[1], SQtrat=f[2],verbose=T)
      AnalAB.=""
      }
    print(AnalAB.)}
  if(verbose==TRUE){ cat("\n")}
}
}


#Desdobramento B/AC
glem=Satterthwaite(QMA = qmea,QMB = qmec,DFa = glea,DFb = glec,b = c)
qmem=(qmea+(c-1)*qmec)/c

glem2=Satterthwaite(QMA = qmeb,QMB = glem,DFa = gleb,DFb = glem,b =a*c )
qmem2=(qmea+(a*c-1)*qmem)/(a*c)


for(i in unique(FA)){
  for(j in unique(FC)){
    id=(FA==i)&(FC==j)







    if(verbose==TRUE){ cat(".................................................................", "\n")}
    if(verbose==TRUE){ cat(paste("Desdobramento do fator B dentro do nivel",i,"do fator A e",j,"do fator C"), "\n")}
    if(verbose==TRUE){ cat(paste("Splitting of factor B within level",i,"of factor A and",j, "of factor C"), "\n")}
    if(verbose==TRUE){
    if(quali[2]==TRUE){
      AnalAB.=ComparacaoMedias(Y[id],FB[id],DFerror = glem2,MSerror = qmem2,alpha = sig)
    }   }

        if(verbose==TRUE){
      if(quali[2]==FALSE){
        f=sq(FB[id],Y[id])
        AnalAB.=RegressaoPolinomial(Y[id],FB[id], glres=glem2, SQres=qmem2*glem2, gltrat=f[1], SQtrat=f[2],verbose=T)
        AnalAB.=""
        }

      print(AnalAB.)}
    if(verbose==TRUE){ cat("\n")}
  }
}


#Desdobramento C/AB
glem=Satterthwaite(QMA = qmea,QMB = qmeb,DFa = glea,DFb = gleb,b = b)
qmem=(qmea+(b-1)*qmeb)/b

glem2=Satterthwaite(QMA = qmec,QMB = glem,DFa = glec,DFb = glem,b =a*b )
qmem2=(qmea+(a*b-1)*qmem)/(a*b)


for(i in unique(FA)){
  for(j in unique(FB)){
    id=(FA==i)&(FB==j)




    if(quali[3]==TRUE){
      AnalAB.=ComparacaoMedias(Y[id],FC[id],DFerror = glem2,MSerror = qmem2,alpha = sig)
    }



    if(verbose==TRUE){ cat(".................................................................", "\n")}
    if(verbose==TRUE){ cat(paste("Desdobramento do fator C dentro do nivel",i,"do fator A e",j,"do fator B"), "\n")}
    if(verbose==TRUE){ cat(paste("Splitting of factor C within level",i,"of factor A and",j, "of factor B"), "\n")}

        if(verbose==TRUE){
      if(quali[3]==FALSE){
        f=sq(FC[id],Y[id])
        AnalAB.=RegressaoPolinomial(Y[id],FC[id], glres=glem2, SQres=qmem2*glem2, gltrat=f[1], SQtrat=f[2],verbose=T)
        AnalAB.=""
        }

      print(AnalAB.)}
    if(verbose==TRUE){ cat("\n")}
  }
}
}









