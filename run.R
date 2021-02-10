rm(list=ls())
library(extraDistr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("funzioni definitive.R")

help(xtable)
#source("C:\\Users\\Baccini\\Documents\\COVID19\\simulazioni_test_rotazione\\funzioni definitive.R")

#screening_alternati(N=N, sens=sens, tempo.sim = 100, R0 = r0, T.exit = t)


################################################# NUOVO CICLO MICHELA CON CALCOLO A 7, 14, 21 E 28 PER TUTTE LE STRATEGIE 
N=24
r0=c( 1.1,2,  3, 5)
t=c(7,14,21) 
comb<-expand.grid(r0,t)
test_A1<-list()
test_A2<-list()
test_B1<-list()
test_B2<-list()
test_C<-list()
test_D<-list()

t1<-Sys.time()
QUANTE.SIM=1000
#s3<-list()
l=1
i=1
p=0.7
for(j in 1: nrow(comb)){
  
  list_lag_A1<- list_lag_A2<- list_lag_B1<- list_lag_B2<- list_lag_C<-list_lag_D<-list()
  for(l in 1:7){
    list_lag_A1[[l]]<-screening_A1(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    list_lag_B1[[l]]<-screening_B1(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    list_lag_D[[l]]<-screening_D(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    print("lag", l)
  }
  for(l in 1:10){
    list_lag_C[[l]]<-screening_C(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    print("lag", l)
  }
  
  for(l in 1:14){
    list_lag_A2[[l]]<-screening_A2(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    list_lag_B2[[l]]<-screening_B2(N=N, sens=p, tempo.sim = 100, R0 = comb[j,1], T.exit = comb[j,2], Nsim=QUANTE.SIM, binom=1, lag=(l-1))
    print("lag", l)
  }
  
  
  test_A1[[j]]<- list_lag_A1  
  test_A2[[j]]<- list_lag_A2
  test_B1[[j]]<- list_lag_B1  
  test_B2[[j]]<- list_lag_B2  
  test_C[[j]] <- list_lag_C
  test_D[[j]] <- list_lag_D
  print(j)
} 

PROB_0pos<-PROB_Q0_10<-PROB_Q0_7<-PROB_Q7_14<-PROB_Q14_21<-PROB_Q10_20<-PROB_Q21_28<-matrix(0, nrow=12, 6)
colnames(PROB_0pos)<-colnames(PROB_Q0_7)<-colnames(PROB_Q0_10)<-colnames(PROB_Q10_20)<-colnames(PROB_Q7_14)<-colnames(PROB_Q14_21)<-colnames(PROB_Q21_28)<-c("A1", "A2", "B1", "B2", "C", "D")

for(j in 1:12){
  PROB_Q0_7[j,]<-c(mean(unlist(lapply(test_A1[[j]], function(x) x$prob_q7))),
                   mean(unlist(lapply(test_A2[[j]],function(x) x$prob_q14))[8:14])*7/14,
                   mean(unlist(lapply(test_B1[[j]], function(x) x$prob_q7))),
                   mean(unlist(lapply(test_B2[[j]],function(x) x$prob_q14))[8:14])*7/14,
                   sum(unlist(lapply(test_C[[j]],function(x) x$prob_q10))[4:10])/10,
                   mean(unlist(lapply(test_D[[j]], function(x) x$prob_q7)))
  )
  PROB_Q7_14[j,]<-c(mean(unlist(lapply(test_A1[[j]], function(x) x$prob_q14))),
                    mean(unlist(lapply(test_A2[[j]], function(x) x$prob_q14)))-mean(unlist(lapply(test_A2[[j]],function(x) x$prob_q14))[8:14])*7/14,
                    mean(unlist(lapply(test_B1[[j]], function(x) x$prob_q14))),
                    mean(unlist(lapply(test_B2[[j]], function(x) x$prob_q14)))-mean(unlist(lapply(test_B2[[j]],function(x) x$prob_q14))[8:14])*7/14,
                    (sum(unlist(lapply(test_C[[j]],function(x) x$prob_q10))[1:3])+                
                       sum(unlist(lapply(test_C[[j]],function(x) x$prob_q20))[7:10]))/10,
                    #NA,
                    mean(unlist(lapply(test_D[[j]], function(x) x$prob_q14)))
  )
  
  PROB_Q14_21[j,]<-c(mean(unlist(lapply(test_A1[[j]], function(x) x$prob_q21))),#vederlo tra 14 e 20 giorni dal primo infetto
                     mean(unlist(lapply(test_A2[[j]],function(x) x$prob_q28))[8:14])*7/14,
                     #NA,
                     mean(unlist(lapply(test_B1[[j]], function(x) x$prob_q21))),
                     #NA,
                     #NA,
                     mean(unlist(lapply(test_B2[[j]],function(x) x$prob_q28))[8:14])*7/14,
                     (sum(unlist(lapply(test_C[[j]],function(x) x$prob_q20))[1:6])+
                        sum(unlist(lapply(test_C[[j]],function(x) x$prob_q30))[10]))/10,
                     mean(unlist(lapply(test_D[[j]], function(x) x$prob_q21)))
  ) 
  
  PROB_Q21_28[j,]<-c(mean(unlist(lapply(test_A1[[j]], function(x) x$prob_q28))),
                     mean(unlist(lapply(test_A2[[j]], function(x) x$prob_q28)))-mean(unlist(lapply(test_A2[[j]],function(x) x$prob_q28))[8:14])*7/14,
                     mean(unlist(lapply(test_B1[[j]], function(x) x$prob_q28))),
                     mean(unlist(lapply(test_B2[[j]], function(x) x$prob_q28)))-mean(unlist(lapply(test_B2[[j]],function(x) x$prob_q28))[8:14])*7/14,
                     #NA
                     sum(unlist(lapply(test_C[[j]],function(x) x$prob_q30))[3:9])/10,
                     mean(unlist(lapply(test_D[[j]], function(x) x$prob_q28)))
  )
  
 # PROB_0pos[j,]<-c(mean(unlist(lapply(test_A1[[j]], function(x) x$prob_0pos))),
 #                    mean(unlist(lapply(test_A2[[j]], function(x) x$prob_0pos))),
  #                  mean(unlist(lapply(test_B1[[j]], function(x) x$prob_0pos))),
   #                  mean(unlist(lapply(test_B2[[j]], function(x) x$prob_0pos))),
   #                  sum(unlist(lapply(test_C[[j]],function(x) x$prob_0pos))),
    #                 mean(unlist(lapply(test_D[[j]], function(x) x$prob_noseen)))
 # )
  
}

#################### se si calcola 1-PROB_Q0_7-PROB_Q7_14-PROB_Q14_21-PROB_Q21_28 SI OTTINE LA PROB DI NESSUN POSITIVO A 28 GG
#################### CHE COINCIDE CON PROB_0pos TRANNE CHE PER LA COLONNA C

#################### GRAFICO MICHELA
load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.7_MB_1000.RData")

#load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws0.9_1000.RData")

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
# j<-4
# plot(MAT.LIST[[j]][1,],ylim=c(0,1),type="l")
# for(k in 2:6){
#   lines(MAT.LIST[[j]][k,],col=k)
# }
# legend(3,0.4,legend=c("A1","A2","B1","B2","C","D"),lty=1,col=1:6)

# 
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
# legend('bottom',legend = c("A1", "A2", "B1", "B2", "C", "D"), col = 1:6, lwd = 3, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
# 
# ########################################

################### da qui non ho rifatto nulla



 

#load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.7_10000.RData")

#load("~/Dropbox/EPM2019/Firenze/COVID Cecilia Giulia/screening scuole/ws_280.9_MB_1000.RData")



c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

#colori<-c25[c(25,1,5,10,7,6)]
colori<-1:6
pie(rep(1,6), col=colori)
 
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


# pdf("densities.pdf")
# par(mfrow=c(3,4))
#  for(i in 1:12){
#     plot(density(infetti_A1[i,]),main=paste("R0=",toString(comb[i,1]),", T=", toString(comb[i,2]) ), ylab="", xlab="")
#    lines(density(infetti_A2[i,]), col=2)
#    lines(density(infetti_B1[i,]), col=3)
#    lines(density(infetti_B2[i,]), col=4)
#    
#    lines(density(infetti_C[i,]), col=5)
#    
#    lines(density(infetti_D[i,]), col=6)
#    legend("topright", legend=c("A1", "A2", "B1", "B2", "C", "D"), lty=1, col=1:6)
#  }
# dev.off()

nome_b<-paste("boxplot_mb_28",  toString(p),"_",toString(QUANTE.SIM),".pdf", sep="" )

pdf(nome_b, width=9, height=10)
par(mfrow=c(3,4))

for(i in ordine_k){
  boxplot(infetti_A1[i,],infetti_A2[i,], infetti_B1[i,], infetti_B2[i,],infetti_C[i,],infetti_D[i,],main=paste("R0=",toString(comb[i,1])," T=", toString(comb[i,2]), "beta=", toString(round(betas[i],2) )), ylab="", xlab="", names=c("A1", "A2", "B1", "B2", "C", "D"), horizontal = TRUE, outline =FALSE, col=c("red","yellow", "yellow" , "green", "green", "green") )
  # legend("topright", legend=c("A1", "A2", "B1", "B2", "C", "D"), lty=1, col=1:6)
}
dev.off()

nome_ws<-paste("ws_28",  toString(p),"_MB_",toString(QUANTE.SIM),".RData", sep="" )

#save.image(nome_ws)







apply(mat_test1, 1,boxplot)

par(mfrow=c(1,1))
boxplot(boxplot_test1, boxplot_test1_pool, boxplot_test2, boxplot_test2_pool,boxplot_test3,names =c("4g, p=0.7","4gpool, p=0.9", "2g, p=0.7" ,"2gpool, p=0.9", "3g, p=0.7"))
print(xtable(test1_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
print(xtable(test1_pool_m, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
print(xtable(test2_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
print(xtable(test2_m_pool, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
print(xtable(test3_m, caption="Screening a 3 gruppi"), include.rownames=FALSE) 



PROB_0pos

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

# 
# matrix_ex_08_3g<-t(matrix(unlist(test1), ncol=nrow(comb)))
# matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:4], 1, cumsum))
# test1_m<-cbind(comb,matr_cumsum_08_3g ,matrix_ex_08_3g[,5:7])
# 
# 
# matrix_ex_08_3g<-t(matrix(unlist(test1_pool), ncol=nrow(comb)))
# matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:4], 1, cumsum))
# test1_pool_m<-cbind(comb,matr_cumsum_08_3g ,matrix_ex_08_3g[,5:7])
# 
# 
# matrix_ex_08_3g<-t(matrix(unlist(test2_pool), ncol=nrow(comb)))
# matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:4], 1, cumsum))
# test2_m_pool<-cbind(comb,matr_cumsum_08_3g ,matrix_ex_08_3g[,5:7])
# 
# matrix_ex_08_3g<-t(matrix(unlist(test2), ncol=nrow(comb)))
# matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:4], 1, cumsum))
# test2_m<-cbind(comb,matr_cumsum_08_3g ,matrix_ex_08_3g[,5:7])
# 
# 
# colnames(test2_m_pool)<-colnames(test2_m)<-colnames(test1_pool_m)<-colnames(test1_m)<-c("R0", "T", "p7", "p14", "p21", "p28", "p0", "media", "varianza")
#  
# 
# matrix_ex_08_3g<-t(matrix(unlist(test3), ncol=nrow(comb)))
# matr_cumsum_08_3g<-t(apply(matrix_ex_08_3g[,1:3], 1, cumsum))
# test3_m<-cbind(comb,matr_cumsum_08_3g ,matrix_ex_08_3g[,4:6])
# colnames(test3_m)<-c("R0", "T", "p10", "p20", "p30",  "p0", "media", "varianza")
# 
# 
# 
# print(xtable(test1_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
# print(xtable(test1_pool_m, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
# print(xtable(test2_m, caption="Screening a 4 gruppi, p=0.7"), include.rownames=FALSE)
# print(xtable(test2_m_pool, caption="Screening a 4 gruppi,con pooltesting, p=0.9"), include.rownames=FALSE)
# print(xtable(test3_m, caption="Screening a 3 gruppi"), include.rownames=FALSE) 
# 
# 
# 
# 
# m_07_3g
# m_08_3g
# m_09_3g
# 
# m_07_p_3g
# m_08_p_3g
# m_09_p_3g
# 
# 
# screening_alternati_ripetuti
# 
# 
# for(i in 1: nrow(comb)){
#   #s1[[i]]<- screening_alternati(N=N, sens=sens, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000)
#   print("c")
#   # s_08_3g[[i]]<- screening_alternati_3g(N=N, sens=0.8, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   # s_09_3g[[i]]<- screening_alternati_3g(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   s_07_3g[[i]]<- screening_alternati_ripetuti(N=N, sens=0.7, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   print("c")# s_09[[i]]<- screening_alternati_3g(N=N, sens=sens, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000)
#   #  s_08_p_3g[[i]]<- screening_alternati_3g(N=N, sens=0.8, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=0)
#   #  s_09_p_3g[[i]]<- screening_alternati_3g(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=0)
#   s_07_p_3g[[i]]<- screening_alternati_ripetuti_pooltesting(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   print(i)
#   
# } 
# 
# 
# matrix_ex_07_3g<-t(matrix(unlist(s_07_3g), ncol=nrow(comb)))
# matr_cumsum_07_3g<-t(apply(matrix_ex_07_3g[,1:4], 1, cumsum))
# m_07_3g<-cbind(comb,matr_cumsum_07_3g ,matrix_ex_07_3g[,5:7])
# colnames(m_07_3g)<-c("R0", "T", "p7", "p14", "p21", "p28", "p0", "media", "varianza")
# 
# 
# 
# matrix_ex_07_pool<-t(matrix(unlist(s_07_p_3g), ncol=nrow(comb)))
# matr_cumsum_07_pool<-t(apply(matrix_ex_07_pool[,1:4], 1, cumsum))
# m_07_3g_pool<-cbind(comb,matr_cumsum_07_pool ,matrix_ex_07_pool[,5:7])
# colnames(m_07_3g_pool)<-c("R0", "T", "p7", "p14", "p21", "p28", "p0", "media", "varianza")
# 
# 
# m_07_p_3g
# m_08_p_3g
# m_09_p_3g
# 
# s_09_pool4g<-list()
#  
# 
# for(i in 1: nrow(comb)){
#   #s1[[i]]<- screening_alternati(N=N, sens=sens, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000)
#   print("c")
#   # s_08_3g[[i]]<- screening_alternati_3g(N=N, sens=0.8, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   # s_09_3g[[i]]<- screening_alternati_3g(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   s_09_pool4g[[i]]<- screening_alternati_4g_pooltesting(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=1)
#   print("c")# s_09[[i]]<- screening_alternati_3g(N=N, sens=sens, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000)
#   #  s_08_p_3g[[i]]<- screening_alternati_3g(N=N, sens=0.8, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=0)
#   #  s_09_p_3g[[i]]<- screening_alternati_3g(N=N, sens=0.9, tempo.sim = 100, R0 = comb[i,1], T.exit = comb[i,2], Nsim=10000, binom=0)
#  
#   
# } 
# 
# 
# 
# 
# 
# matrix_ex_09_pool_4g<-t(matrix(unlist(s_09_pool4g), ncol=nrow(comb)))
# matr_cumsum_09_pool_4g<-t(apply(matrix_ex_09_pool_4g[,1:4], 1, cumsum))
# m_07_09_pool_4g<-cbind(comb,matr_cumsum_09_pool_4g ,matrix_ex_09_pool_4g[,5:7])
# colnames(m_07_09_pool_4g)<-c("R0", "T", "p7", "p14", "p21", "p28", "p0", "media", "varianza")
# 
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

  
#par(mfrow=c(4,4))

#n <- 6
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

