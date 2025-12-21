kalman.miss.wge=function(y,start,gam0,F,gamV,Gtmiss,gamW){
  n=length(y)
  sqgamV=sqrt(gamV)
  sqgamW=sqrt(gamW)
  chgamW=chol(gamW)
  yy=y
  #Fix from David Stoffer 1-29-2023
  result=astsa::Ksmooth(y=yy, A=Gtmiss, mu0=start, Sigma0=gam0, Phi=F, sQ=sqgamV, sR=chgamW)
  # end fix
  
  #Bivin fixed xp, xf and xs to Xp, Xf and Xs below 1-20-2023
  pfs=cbind(yy,result$Xp,result$Pp,result$Xf,result$Pf,result$Xs,result$Ps)
  colnames(pfs)=c("Data","Prediction","Var_Predict","Filter","Var_Filter","Smooth","Var_Smooth")
  return(pfs)
}