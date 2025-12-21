kalman.wge=function(y,start,gam0,F,gamV,G,gamW){
  n=length(y)
  sqgamV=sqrt(gamV)
  sqgamW=sqrt(gamW)
  yy=y
  #Fix from David Stoffer 1-29-2023
  result=astsa::Ksmooth(y=yy, A=G, mu0=start, Sigma0=gam0, Phi=F, sQ=sqgamV, sR=sqgamW)
  pfs=cbind(yy,result$Xp,result$Pp,result$Xf,result$Pf,result$Xs,result$Ps)
  colnames(pfs)=c("Data","Prediction","Var_Predict","Filter","Var_Filter","Smooth","Var_Smooth")
  return(pfs)
}
