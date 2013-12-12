library(MASS); library(mvtnorm); library(survival)

logLam0.true = function(tt){2*log(tt)-3};

gam.true = t(c(0.5))
yes.DepCen = 1 < 0; 
yes.Z = 1<0
Zmatch = 1<0
N = 2000; 
B0 = 200
Zmatch.ind = NULL 
t0 =1; 
type ='cutoff'
uu0 = c(0        )
m.match = 3
myind.Zmatch = NULL

mydata.cht = SIM.FUN(nn=N,DepCen=yes.DepCen,Zmatch=yes.Z,a0=a0.match)
mydata.ncc = mydata.cht[[1]];

Vij.IND = mydata.cht[[2]]; 
Iik0.mat = mydata.cht[[3]]

junk=NCC.Cox.explicit.PTB.s.y(data=mydata.ncc,V.IND=Vij.IND,Iik0=Iik0.mat,wgtk.ptb=NULL,B0=500,Zmatch.ind=myind.Zmatch, t0 = 1);

Ptb.result = matrix(0,B0,1+4*length(uu0))

for (i in 1:B0) {Ptb.result[i,]=unlist(Ptb.ROC.FUN(junk$ck,
                                                   junk$CondSck.ptb[,i],
                                                   junk$Fck.ptb[,i],uu0,type, 
                                                   junk$yk))
}


c(junk$beta.est, apply(Ptb.result,2,mean))
c(junk$beta.sd,apply(Ptb.result,2,sd))



mydata.cht <- SIM.data.singleMarker( nn=2000, 
                                     beta = log(3), 
                                     lam0 = .1, 
                                     cens.perc = .8, 
                                     time.max = NULL)

getEstandSE_SP(data = mydata.cht, 
               cutpoint = 0,
               predict.time = 1)
