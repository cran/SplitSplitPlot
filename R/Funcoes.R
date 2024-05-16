###################################################
Satterthwaite=function(QMA,QMB,DFa,DFb,b) {(QMA+(b-1)*QMB)^2/((QMA^2/DFa)+(((b-1)*QMB)^2)/DFb)}
##################################################################################


sq=function(Fator,Y){
  X=as.factor(Fator)
  #X=c(X)
  NumTrat=length(unique(as.factor(Fator)));
  NomeTrat=unique(as.factor(Fator))
  NumParc=length(X);

  sq=sum(tapply(Y,X,sum)^2/table(X))-sum(Y)^2/NumParc


  return(c(GL=NumTrat-1,SQ=sq))
}
##############################################################
sqI2=function(Fator1,Fator2,Y){
  Fator1=as.factor(Fator1)
  Fator2=as.factor(Fator2)
  X=paste(Fator1,Fator2)
  sqA=sq(Fator1,Y)
  sqB=sq(Fator2,Y)
  sqA.B=sq(X,Y)

  return(sqA.B-sqA-sqB)
}
##################################################################


##############################################################
sqI3=function(Fator1,Fator2,Fator3,Y){
  Fator1=as.factor(Fator1)
  Fator2=as.factor(Fator2)
  Fator3=as.factor(Fator3)
  X=paste(Fator1,Fator2,Fator3)
  sqA=sq(Fator1,Y)
  sqB=sq(Fator2,Y)
  sqC=sq(Fator3,Y)

  sqAB=sqI2(Fator1,Fator2,Y)
  sqAC=sqI2(Fator1,Fator3,Y)
  sqBC=sqI2(Fator2,Fator3,Y)

  sqA.B.C=sq(X,Y)

  return(sqA.B.C-sqA-sqB-sqC-sqAB-sqAC-sqBC)
}
##################################################################
