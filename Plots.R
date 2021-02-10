


####################
load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.7_MB_1000.RData")

 
nome_p<-paste("probas_",  toString(p),"_",toString(QUANTE.SIM),".pdf", sep="" )
pdf(nome_p, height=12,width=9)
n <- 6
par(omi=c(0.4,0.4,0.4,0.4),mfrow = c(3, 4), mar=c(4, 4, 3,1.5) )
betas<-comb[,1]/comb[,2]
MAT.LIST<-list()
for(j in 1:12){
  yax=" "
  if(j %in% c(1,5,9)){
    yax="p"
  }
  
  MATRICE<-NULL
  for (k in 1:6){
    vec<-c(PROB_Q0_7[j,k],PROB_Q7_14[j,k],PROB_Q14_21[j,k],PROB_Q21_28[j,k])
    MATRICE<-rbind(MATRICE,cumsum(vec))}
  MAT.LIST[[j]]<-MATRICE}

par(mfrow=c(3,4))

ordine_k<-c(9,10,11,12,5,6,7,8,1,2,3,4)

for(j in ordine_k){
  plot(c(7,14,21,28),MAT.LIST[[j]][1,],ylim=c(0,1),xlab="t",xaxt="n",lwd=2,type="b",main=paste("R0=",toString(comb[j,1])," T=", toString(comb[j,2]), "beta=", toString(round(betas[j],2) )), ylab=yax, col=1)
  axis(1, at=c(7,14,21,28), labels=c(7,14,21,28))
  for(k in 2:6){
    lines(c(7,14,21,28),MAT.LIST[[j]][k,],col=k,type="b",lwd=2)
  }
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("A1", "A2", "B1", "B2", "C", "D"), col = 1:6, lwd = 3, xpd = TRUE, horiz = TRUE, cex = 2, seg.len=1, bty = 'n')

dev.off()
 

 

 
 

infetti_A1<-infetti_B1<-infetti_D<-matrix(0,nrow(comb) ,QUANTE.SIM*7)
infetti_A2<-infetti_B2<-matrix(0,nrow(comb) ,QUANTE.SIM*14)
infetti_C<-matrix(0,nrow(comb) ,QUANTE.SIM*10)

for(j in 1:nrow(comb)){
  infetti_A1[j,]<-unlist(lapply(test_A1[[j]], function(x) x$infetti.fino.a.qui))
  infetti_A2[j,]<-unlist(lapply(test_A2[[j]], function(x) x$infetti.fino.a.qui))
  infetti_B1[j,]<-unlist(lapply(test_B1[[j]], function(x) x$infetti.fino.a.qui))
  infetti_B2[j,]<-unlist(lapply(test_B2[[j]], function(x) x$infetti.fino.a.qui))
  infetti_C[j,]<-unlist(lapply(test_C[[j]], function(x) x$infetti.fino.a.qui))
  infetti_D[j,]<-unlist(lapply(test_D[[j]], function(x) x$infetti.fino.a.qui))
}

 
nome_b<-paste("boxplot_mb_28",  toString(p),"_",toString(QUANTE.SIM),".pdf", sep="" )

pdf(nome_b, width=9, height=10)
par(mfrow=c(3,4))

for(i in ordine_k){
  boxplot(infetti_A1[i,],infetti_A2[i,], infetti_B1[i,], infetti_B2[i,],infetti_C[i,],infetti_D[i,],main=paste("R0=",toString(comb[i,1])," T=", toString(comb[i,2]), "beta=", toString(round(betas[i],2) )), ylab="", xlab="", names=c("A1", "A2", "B1", "B2", "C", "D"), horizontal = TRUE, outline =FALSE, col=c("red","yellow", "yellow" , "green", "green", "green") )
  }
dev.off()

nome_ws<-paste("ws_28",  toString(p),"_MB_",toString(QUANTE.SIM),".RData", sep="" )

  

apply(mat_test1, 1,boxplot)

par(mfrow=c(1,1))
boxplot(boxplot_test1, boxplot_test1_pool, boxplot_test2, boxplot_test2_pool,boxplot_test3,names =c("4g, p=0.7","4gpool, p=0.9", "2g, p=0.7" ,"2gpool, p=0.9", "3g, p=0.7"))
print(xtable(test1_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
print(xtable(test1_pool_m, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
print(xtable(test2_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
print(xtable(test2_m_pool, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
print(xtable(test3_m, caption="Screening a 3 gruppi"), include.rownames=FALSE) 



 

############tabella 1#####################

PROB_0pos<-1-PROB_Q0_7-PROB_Q7_14-PROB_Q14_21-PROB_Q21_28 

mat<-cbind(comb,betas,PROB_0pos)
colnames(mat)<-c("R0", "T","beta","A1", "A2", "B1", "B2", "C", "D")
mat0<-mat[1:4,]
mat1<-mat[9:12,]
mat[1:4,]<-mat1
mat[9:12,]<-mat0
print(xtable(mat, caption="", digits = 3), include.rownames=FALSE) 

#################tabella 2####################3

mat2<-cbind(comb,betas,apply(infetti_A1, 1, mean),
            apply(infetti_A2, 1, mean),
            apply(infetti_B1, 1, mean),
            apply(infetti_B2, 1, mean),
            apply(infetti_C, 1, mean),
            apply(infetti_D, 1, mean))
colnames(mat2)<-c("R0", "T","beta","A1", "A2", "B1", "B2", "C", "D")
mat0<-mat2[1:4,]
mat1<-mat2[9:12,]
mat2[1:4,]<-mat1
mat2[9:12,]<-mat0

print(xtable(mat2, caption="Screening a 3 gruppi"), include.rownames=FALSE) 

fq<-function(x){
  return(quantile(x,0.90))
}


mat3<-cbind(comb,betas,apply(infetti_A1, 1, fq),
            apply(infetti_A2, 1, fq),
            apply(infetti_B1, 1, fq),
            apply(infetti_B2, 1, fq),
            apply(infetti_C, 1, fq),
            apply(infetti_D, 1, fq))
colnames(mat3)<-c("R0", "T","beta","A1", "A2", "B1", "B2", "C", "D")
mat0<-mat3[1:4,]
mat1<-mat3[9:12,]
mat3[1:4,]<-mat1
mat3[9:12,]<-mat0
print(xtable(mat, caption="Probabity of not seeing the infection after 35 days", digits=c(1, 1,0, 2,rep(3,6))), include.rownames=FALSE) 
print(xtable(mat2, caption="Mean of the lost infection days",digits=c(1, 1,0,2, rep(1,6))), include.rownames=FALSE) 
print(xtable(mat3, caption="90% quantile of the lost infecton days",digits=c(1, 1,0,2, rep(0,6))), include.rownames=FALSE) 


fq<-function(x){
  return(quantile(x,0.90))
}


mat4<-cbind(comb,betas
            ,apply(infetti_A1, 1, mean),apply(infetti_A1, 1, fq),
            apply(infetti_A2, 1, mean),apply(infetti_A2, 1, fq),
            apply(infetti_B1, 1, mean),apply(infetti_B1, 1, fq),
            apply(infetti_B2, 1, mean),apply(infetti_B2, 1, fq),
            apply(infetti_C, 1, mean),apply(infetti_C, 1, fq),
            apply(infetti_D, 1, mean),apply(infetti_D, 1, fq))
#colnames(mat3)<-c("R0", "T","A1", "A2", "B1", "B2", "C", "D")

mat0<-mat4[1:4,]
mat1<-mat4[9:12,]
mat4[1:4,]<-mat1
mat4[9:12,]<-mat0
colnames(mat4)<-c(" ", " "," ",rep(c("mean", "90%"),6 ))
print(xtable(mat4, caption="", digits=c(1, 1,0, 2,rep(c(1,0),6))), include.rownames=FALSE) 
  matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:3], 1, cumsum))
 
 
load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.7_MB_1000.RData")


MAT.LIST<-list()
for(j in 1:12){
  MATRICE<-NULL
  for (k in 1:6){
    vec<-c(PROB_Q0_7[j,k],PROB_Q7_14[j,k],PROB_Q14_21[j,k],PROB_Q21_28[j,k])
    MATRICE<-rbind(MATRICE,cumsum(vec))}
  MAT.LIST[[j]]<-MATRICE}

par(mfrow=c(3,4))
MATLIST07<-MAT.LIST

load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.9_MB_1000.RData")

MAT.LIST<-list()
for(j in 1:12){
  MATRICE<-NULL
  for (k in 1:6){
    vec<-c(PROB_Q0_7[j,k],PROB_Q7_14[j,k],PROB_Q14_21[j,k],PROB_Q21_28[j,k])
    MATRICE<-rbind(MATRICE,cumsum(vec))}
  MAT.LIST[[j]]<-MATRICE}

par(mfrow=c(3,4))

MATLIST09<-MAT.LIST


 
betas<-comb[,1]/comb[,2]


ordine_k2<-c(9,4)
nome_p<-paste("probas_2",  toString(p),"_",toString(QUANTE.SIM),".pdf", sep="" )
pdf(nome_p, height=11,width=16)


par( mfrow = c(1, 2), omi=c(1.5,0,0,0))

for(j in ordine_k2){
  plot(c(7,14,21,28),MATLIST07[[j]][1,],ylim=c(0,1),xaxt="n", lwd=2,type="b",cex.main = 2, main=paste("R0=",toString(comb[j,1])," T=", toString(comb[j,2]), "beta=", toString(round(betas[j],2) )), xlab="t" ,ylab="", col=1)
  axis(1, at=c(7,14,21,28), labels=c(7,14,21,28))
  
  lines(c(7,14,21,28),MATLIST09[[j]][1,],col=1,type="b",lwd=2, lty=2)
  # legend('bottom',legend = c("A1 p=0.7", "B1 p=0.7", "D p=0.7", "A1 p=0.9", "B1 p=0.9","D p=0.9"), col = c(1,3,6,1,3,6), lwd = 3,lty=c(1,1,1,2,2,2), xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
  
  for(k in c(3,6)){
    lines(c(7,14,21,28),MATLIST07[[j]][k,],col=k,type="b",lwd=2)
    lines(c(7,14,21,28),MATLIST09[[j]][k,],col=k,type="b",lwd=2, lty=2)
  }
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 3, 0), mar = c(0, 0, 5, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("A1 p=0.7", "B1 p=0.7", "D p=0.7", "A1 p=0.9", "B1 p=0.9","D p=0.9"), col = c(1,3,6,1,3,6), lwd = 2.3,lty=c(1,1,1,2,2,2), xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1.5, bty = 'n')

dev.off()