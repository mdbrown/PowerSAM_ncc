
	## Iik0: matrix with nrow=sum(di),ncol=nn-sum(di) ## indicates whether xk is in the matched risk set for ti
	
### for NCC study requires two data set: 
### data.cht: cohort data; with xi,di,vi,yi,zi,id;
### data.matchinfo: caseid, matched controlid; 
	
##Vij.IND: 
	## tivec: failure time of selected cases; 
	## nivec: size of riskset for the case sum(xi>=ti)
	## id of the case; 
	## ids of the matched controls for the case; 
	
NCC.main.FUN = function(data.cht,data.matchinfo,Zmatch,Zmatch.ind=myind.Zmatch,m.match,B0,t0,uu0,type) {
	#mydata.cht = prepareNCCdata(data.cht,data.matchinfo,Zmatch,m.match)
	mydata.cht = SIM.FUN(nn=N,DepCen=yes.DepCen,Zmatch=yes.Z,a0=a0.match)
	mydata.ncc = mydata.cht[[1]];
	Vij.IND = mydata.cht[[2]]; Iik0.mat = mydata.cht[[3]]
	junk=NCC.NP.kn.PTB.s.y(data=mydata.ncc,V.IND=Vij.IND,Iik0=Iik0.mat,wgtk.ptb=NULL,B0=500,Zmatch.ind=myind.Zmatch, t0 = 1);
	
	Ptb.result = matrix(0,B0,1+3*length(uu0))
	
	for (b in 1:B0) {Ptb.result[b,]=unlist(Ptb.ROC.FUN(junk$ck,junk$CondSck.ptb[,b],junk$Fck.ptb[,b],uu0,type))
	}
	keep = apply(!is.na(Ptb.result),1,sum)==dim(Ptb.result)[2]
	cbind(c(junk$beta.est, unlist(Ptb.ROC.FUN(junk$ck,junk$CondSck,junk$Fck,uu0,type))),c(junk$beta.sd,apply(Ptb.result[keep,],2,sd)))
}
	
## data: xi, di, vi, yi; 
NCC.NP.kn.PTB.s.y <- function(data,V.IND,Iik0=Iik0.mat,wgtk.ptb=NULL,B0=500,Zmatch.ind,t0, m.match)
  {     
    bw.power = 0.3
    xi = data[,1]; di = data[,2]; vi = data[,3]; yi = data[,-(1:3),drop=F];yi=as.matrix(yi); n.t0 = length(t0); nv = sum(vi)
    NN = length(xi); nv = sum(vi); tj = V.IND[,1]; pi.tj = V.IND[,2]/NN
    wgtk = 1/FNCC.WGT.FUN(data,V.IND,Iik0=Iik0.mat,Zmatch.ind, m.match = m.match); 
    
     xk = xi[vi==1]; dk = di[vi==1]; yk = yi[vi==1,,drop=F]; idk = (1:NN)[vi==1];  
    yk.old <- yk
    #wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,Zmatch.ind=Zmatch.ind,B0=B0); wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
    wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,B0=B0, m.match=m.match); 
    wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
    npy = dim(yk)[2]
    Shat.yk.ptb=matrix(0,nv,B0)
    
    if (npy>1) {
    	betahat = coxph(Surv(xk,dk)~yk,weights=wgtk,robust=T); 
      betahat=betahat$coef; 
      ebyk = c(exp(yk%*%betahat));
      
     	betaptb = PTB.beta.FUN.explicit(wgtk.ptb,cbind(xk,dk,wgtk,yk),betahat=betahat) 
     	yk.ptb= matrix(0,nv,B0)
     	yk.ptb = yk%*%betaptb
     	yk = ebyk   	
    }
     


     hhat.v = 1.06*min(sd(yk),IQR(yk)/1.34)*nv^(-bw.power)
     Shat.yk = CondSurv.FUN(wgtk,xk, dk, yk, t0, hhat.v); 
     if (npy>1) {
     	for (b in 1:B0) {
     		hhat.v = 1.06*min(sd(yk.ptb[,b]),IQR(yk.ptb[,b])/1.34)*nv^(-bw.power)
        	Shat.yk.ptb[,b] = CondSurv.FUN(wgtk.ptb[,b],xk, dk, yk.ptb[,b], t0, hhat.v); }
     } else {
     		for (b in 1:B0) {
        	Shat.yk.ptb[,b] = CondSurv.FUN(wgtk.ptb[,b],xk, dk, yk, t0, hhat.v); }
     }
      
    Fyk = sum.I(yk, ">=", yk, wgtk)/sum(wgtk) 
    Fyk.ptb = matrix(0,length(yk),B0)
  
    for (b in 1:B0) {Fyk.ptb[,b] =sum.I(yk, ">=", yk, wgtk.ptb[,b])/sum(wgtk.ptb[,b])}
    #Fyk.ptb = apply(wgtk.ptb,2,myfyk(yk,wgtk))
    if (npy>1) {
	list(ck =yk, 
       beta.est=betahat,
       beta.sd =apply(betaptb,1,sd), 
       CondSck = Shat.yk,
       CondSck.ptb = Shat.yk.ptb, 
       Fck = Fyk,
       Fck.ptb = Fyk.ptb, 
       yk = yk.old)
    } else{ 
      list(ck =yk, 
           CondSck = Shat.yk,
           CondSck.ptb = Shat.yk.ptb, 
           Fck = Fyk,
           Fck.ptb = Fyk.ptb, 
           yk = yk.old)  }
}
 

 PTB.beta.FUN.explicit <- function(wgt.ptb,data,betahat)
  { 
    xi = data[,1]; di = data[,2]; wgti = data[,3]; yi = data[,-(1:3),drop=F]; ebyi = c(exp(yi%*%betahat)); py = ncol(yi)
    tmpind = di==1; tj = xi[tmpind]; wgt.ptb.j = wgt.ptb[tmpind,]; pv = ncol(wgt.ptb)
    pi0.tj   = c(PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=0))/sum(wgti) 
    pi1.tj   =   PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=1)/sum(wgti)
    Ahat = PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=2)/sum(wgti)*pi0.tj - pi1.tj[,rep(1:py,py)]*pi1.tj[,rep(1:py,rep(py,py))]
    Ahat = matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=py)
    term1 = t(wgt.ptb.j)%*%yi[tmpind,]
    term2 = t(wgt.ptb.j)%*%(pi1.tj/pi0.tj)
    
    pi1.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=1)/sum(wgti) ## ny x nptb
    pi0.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=0)/sum(wgti);## nptb x 1 
    term3 = t(rep(1,length(tj)))%*%((pi1.tj.ptb*pi0.tj-pi0.tj.ptb[,rep(1:pv,rep(py,pv))]*pi1.tj[,rep(1:py,pv)])/pi0.tj^2)
    term3 = t(matrix(term3,nrow=py))
    betastar = betahat + solve(Ahat)%*%t(term1 - term2 - term3)  # n.y x nptb     
    betastar
  }

CondSurv.FUN<-function(IPW, xi, di, yi, tt0, bw)
  {
    nv = length(xi); Mi = rep(0,nv); yy = yi
    kerni.yy = Kern.FUN(c(yy),c(yi),bw)*IPW; ## nv x ny matrix
    skern.yy = colSums(kerni.yy)
    tmpind = (xi<=tt0)&(di==1); tj = xi[tmpind]; nj = length(tj)
    pi.tj.yy = sum.I(tj,"<=",xi,kerni.yy)/VTM(skern.yy,nj) ## nj x ny matrix ##    
    dLam.tj.yy = kerni.yy[tmpind,]/pi.tj.yy/VTM(skern.yy,nj); 
    dLam.tj.yy[is.na(dLam.tj.yy)] = 0
    Shat.t0.yi = exp(-colSums(dLam.tj.yy))
    Shat.t0.yi    
  }
  
  
Kern.FUN <- function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
  { 

    out = (VTM(zz,length(zi))- zi)/bw
    switch(kern0,
            "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
            "gauss"= dnorm(out)/bw)
  }
  
myfyk=function(yk, wgtk) {sum.I(yk, ">=", yk, wgtk)/sum(wgtk)}
