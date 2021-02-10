
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################

SIR<-function(I.0=1,S.0=20, R.0=0,R0=2.5,Texit=21, N=100,binom=0){
  I<-S<-R<-I.new<-R.new<-NA
  alpha<-1/Texit
  beta <-R0/Texit
  I[1]<-I.0
  S[1]<-S.0-1  
  R[1]<-R.0
 
  for(i in 2:N){
    I.new[i]<-min(S[i-1], rpois(1,beta*I[i-1]*max(0,S[i-1])/S[1]))
    R.new[i]<-min(I[i-1],rpois(1,I[i-1]*alpha))
    
    if(binom==1){
      
      p_SI <- 1 - exp(-beta * I[i-1] / S[1]) # S to I
      p_IR <- 1 - exp(-alpha) # I to R
      ## Draws from binomial distributions for numbers changing between
      ## compartments:
      I.new[i] <- rbinom(1,S[i-1], p_SI)
      R.new[i] <- rbinom(1,I[i-1], p_IR)
    }
    I[i]<-I[i-1]+I.new[i]-R.new[i]
    R[i]<-R[i-1]+R.new[i]
    S[i]<-S[i-1]-I.new[i]
  }
  
  
  I.new[1]=I.0
  R.new[1]=R.0
  
  
  # h=length(I)
  #  I=I[-h]
  #   S=S[-h]
  #  R=R[-h]
  #   R.new=R.new[-h]
  #  I.new=I.new[-h]
  #}
  
  return(res=list(I=I,R=R,S=S,I.new=I.new ,R.new=R.new ))
}

#############################################################################################################
 

d<-SIR(R0=1.1, Texit = 14)
d$I+d$R+d$S
d$I.new
d$S

#####################################################################
# scenario 1: R0=3 in 14 giorni
infetti.sir<-SIR(1,20,0,3,14,60)[[1]]

#############################################
### SIMULAZIONE
#############################################
#  
# screening_alternati<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=0, lag=0){
#   k7.vec<-k14.vec<-pos7.vec<-pos14.vec<-NULL
#   infetti.fino.a.qui<-NULL
#   for (i in 1:Nsim){
#     sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom=binom)
#     infetti<-c(rep(0,lag),sir[[1]])
#     
#     Inew<-c(rep(0,lag),sir$I.new)
#     infetti7<-sum(Inew[1:7])  #ho messo uno perche si inizia sempre con un infetto
#     infetti14<-sum(Inew[1:14])
#     infetti_all<-sum(Inew[1:(35+lag)])
#     ##### QUI NON SONO SICURA CHE VADA AGGIUNTO 1
# 
#   
#  #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato
# 
# 
# k7<-rhyper(nn=1,m=infetti7, n=N-infetti7,k=N/2) #infetti nel primo gruppo 
# k14<-infetti7-k7+ rhyper(nn=1,m=infetti14-infetti7, n=N-infetti14 ,k=N/2-infetti7+k7)
# 
# pos7<-rbinom(1,k7,sens)
# pos14<-rbinom(1,k14,sens)
# k7.vec<-c(k7.vec,k7)
# k14.vec<-c(k14.vec,k14)
# pos7.vec<-c(pos7.vec,pos7)
# pos14.vec<-c(pos14.vec,pos14)
# infetti.fino.a.qui[i]<-sum(infetti[1:7])*(pos7!=0)+sum(infetti[1:14])*(pos14!=0&&pos7==0)+sum(infetti[1:(35+lag)])*(pos14==0&&pos7==0)
# 
# }
# 
# tavola<-table(pos7.vec,pos14.vec)
# 
# #IN QUANTE SITUAZIONI ARRIVO AL SECONDO SCREENING
# prob_2s<-sum(tavola[1,-1])/Nsim
# #IN QUANTE SITUAZIONI NON VEDO I POSITIVI
# prob_0pos<-tavola[1,1]/Nsim
# #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
# prob_q7<-sum(tavola[-1,])/Nsim
# 
# #GIORNI INFEZIONE CHE NON VEDO 
# media<-mean(infetti.fino.a.qui)
# varianza<-sd(infetti.fino.a.qui)^2 
# return(list(prob_q7=prob_q7,prob_2s=prob_2s ,prob_0pos=prob_0pos, media=media, varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
# }


##########3##########3##########3##########3##########3##########3
#                 screening alternati a 4 gruppi 
##########3##########3##########3##########3##########3##########3

########### PIANO D

#### MODIFICATA CON LAG
screening_D<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10, binom=0, lag=0 ){

#screening  alternati su 1/4 della classe a 7, 14, 21, 28

    # N=20
    # sens=0.7
    #  R0=3
    # tempo.sim=60
    # T.exit=14
    # Nsim=1000 
    # lag (da 0 a 6) indica che l'epidemia ? partita da 7-lag giorni
  
  
pos7<-pos14<-pos21<-pos28<-NULL
  g1<-g2<-g3<-g4<-matrix(0, Nsim,5)
  #set.seed(123)
  infetti.fino.a.qui<-NULL
  for (i in 1:Nsim){
    sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom)
    infetti<-c(rep(0,lag),sir[[1]])
    Inew<-c(rep(0,lag),sir$I.new)
    i7<-sum(Inew[1:7])  #ho messo uno perche si inizia sempre con un infetto
    #print(i7)
    i14<-sum(Inew[1:14])
    i21<-sum(Inew[1:21])
    i28<-sum(Inew[1:28])
    iall<-sum(Inew[1:(35+lag)])
##### sommo fino a 1 settimana dopo l'ultimo screening
    #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato
    
      g1[i,]<- rmvhyper(1, c(i7, i14-i7, i21-i14, i28-i21, N-i28), k=N/4)
      g2[i,]<- rmvhyper(1, c(i7, i14-i7, i21-i14, i28-i21, N-i28)-g1[i,], k=N/4)
      g3[i,]<- rmvhyper(1, c(i7, i14-i7, i21-i14, i28-i21, N-i28)-g1[i,]-g2[i,], k=N/4)
      g4[i,]<- c(i7, i14-i7,i21-i14,i28-i21, N-i28 )-g1[i,]-g2[i,]-g3[i,]
      
    pos7[i]<-rbinom(1,g1[i,1],sens)
    pos14[i]<-rbinom(1,g2[i,1]+g2[i,2],sens)
    pos21[i]<-rbinom(1,g3[i,1]+g3[i,2]+g3[i,3],sens)
    pos28[i]<-rbinom(1,g4[i,1]+g4[i,2]+g4[i,3]+g4[i,4],sens)
    
    infetti.fino.a.qui[i]<-sum(infetti[1:7])* (pos7[i]!=0)+
    sum(infetti[1:14])*(pos7[i]==0&&pos14[i]!=0)+
    sum(infetti[1:21])*(pos7[i]==0&&pos14[i]==0&&pos21[i]!=0)+
    sum(infetti[1:28])*(pos7[i]==0&&pos14[i]==0&&pos21[i]==0&&pos28[i]!=0)+
    sum(infetti[1:(28+lag)])*(pos7[i]==0&&pos14[i]==0&&pos21[i]==0&&pos28[i]==0)
    
  }
  
   
  
  #prob_0pos<-sum(t0[1,])/Nsim
  
   
  
  #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
  prob_q7<-sum(pos7>0)/Nsim
 
   

index_no7<-which(pos7==0)
index_no14<-which(pos14==0)
index_no21<-which(pos21==0)
index_no28<-which(pos28==0)
index_no7_no14<-intersect(index_no7,index_no14)
index_no7_no14_no21<-intersect(index_no7_no14, index_no21)
index_no7_no14_no21_no28<-intersect(index_no7_no14_no21, index_no28)
#IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
prob_q7<-sum(pos7>0)/Nsim
#IN QUANTI CASI QUARANTENO DOPO 2 SETTIMANA
prob_q14<-sum(pos14[index_no7]>0)/Nsim
#IN QUANTI CASI QUARANTENO DOPO 3 SETTIMANA
prob_q21<-sum(pos21[index_no7_no14]>0)/Nsim
#IN QUANTI CASI QUARANTENO DOPO 4 SETTIMANA
prob_q28<-sum(pos28[index_no7_no14_no21]>0)/Nsim

#IN QUANTI CASI non quaranteno 
prob_noseen<-length(index_no7_no14_no21_no28)/Nsim



prob_q7+prob_q14+prob_q21+prob_q28+prob_noseen
 
    
  #GIORNI INFEZIONE CHE NON VEDO
  
#  media<-sum(infetti[1:14])*prob_q14+ sum(infetti[1:7])*prob_q7+  sum(infetti[1:21])*prob_q21+sum( infetti[1:28])*prob_q28 +sum(infetti)*prob_noseen
  media=mean(infetti.fino.a.qui)
  #varianza<-(sum(infetti[1:14] )-media)^2*(prob_q14)+
  #  (sum( infetti[1:7])-media)^2*prob_q7+
  #  (sum( infetti[1:21])-media)^2*prob_q21+
   # (sum( infetti[1:28])-media)^2*prob_q28+
   # (sum(infetti)-media)^2*prob_noseen
    varianza=sd(infetti.fino.a.qui)^2
  return(list( prob_q7=prob_q7,prob_q14=prob_q14, prob_q21=prob_q21, prob_q28=prob_q28,prob_noseen=prob_noseen, media=media, varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui, lag=lag))
}


###### SCREENING A TAPPETO A 15 GG
 
#   
# 
# screening_a_tappeto<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, lag=0){
# 
# #### screening su tutti a 14 giorni
# 
#   infetti.fino.a.qui<-NULL
#   
#   pos14T.vec<-NULL
#   for(j in 1:Nsim){
#   sir<-SIR(1,N,0,R0,T.exit,tempo.sim)
#   infetti<-c(rep(0,lag),sir[[1]])
#   
#  Inew<-c(rep(0,lag),sir$I.new) 
#   infetti14<-sum(Inew[1:14])+1
#   infetti_all<-sum(Inew[1:(lag+35)])+1
# pos14<-rbinom(1,infetti14,sens)
# pos14T.vec<-c(pos14T.vec,pos14)
# infetti.fino.a.qui[j]<-sum(infetti[1:14])*(pos14!=0)+sum(infetti[1:(35+lag)])*(pos14==0)
# }
# tavolaT<-table(pos14T.vec)
# #IN QUANTE SITUAZIONI NON VEDO I POSITIVI
# prob_0pos<-as.numeric(tavolaT[1]/Nsim)
# 
# #GIORNI INFEZIONE CHE NON VEDO
# #media<-sum(round(infetti[1:14],0))*(1-prob_0pos)+ sum(round(infetti[1:21],0))*prob_0pos
# media<-mean(infetti.fino.a.qui)
# varianza<-sd(infetti.fino.a.qui)^2
# #varianza<-(sum(round(infetti[1:14],0))-media)^2*(1-prob_0pos)+ (sum(round(infetti[1:21],0))-media)^2*(1-prob_0pos)
# return(list(tavola=tavolaT, prob_0pos=prob_0pos,media=media , varianza=varianza))
# }
# 
# 
# screening_a_tappeto()


#############################################
############################################# MICHELA
#############################################

### 3 gruppi, cadenza 10, 20, 30 gg

##########3##########3##########3##########3##########3##########3
#                 screening alternati a 3 gruppi 
##########3##########3##########3##########3##########3##########3
# MODIFICATA CON LAG

#### PIANO C


screening_C<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=0, lag=0){

# screening su 1/3 della classe a 10, 20, 30

  # 
  # N=20
  # sens=0.7
  #  R0=3
  # tempo.sim=60
  # T.exit=14
  # Nsim=1000
  # binom=1
  
  
  pos10<-pos20<-pos30<-NULL
  g1<-g2<-g3<-matrix(0, Nsim,4)
  #set.seed(123)
  infetti.fino.a.qui<-NULL
  for (i in 1:Nsim){
    sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom)
    infetti<-c(rep(0,lag),sir[[1]])
    
Inew<-c(rep(0,lag),sir$I.new)    
i10<-sum(Inew[1:10])  #ho messo uno perche si inizia sempre con un infetto
    i20<-sum(Inew[1:20])
    i30<-sum(Inew[1:30])
    iall<-sum(Inew[1:(lag+35)])
    #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato
    
    N1<-round(N/3)
    N2<-round(N/3)
    N3<- N-N1-N2   
    
    g1[i,]<- rmvhyper(1, c(i10, i20-i10, i30-i20, N-i30), k=N1)
    g2[i,]<- rmvhyper(1, c(i10, i20-i10, i30-i20, N-i30)-g1[i,], k=N2)
    g3[i,]<- c(i10, i20-i10, i30-i20, N-i30)-g1[i,]-g2[i,]
    
    
    pos10[i]<-rbinom(1,g1[i,1],sens)
    pos20[i]<-rbinom(1,g2[i,1]+g2[i,2],sens)
    pos30[i]<-rbinom(1,g3[i,1]+g3[i,2]+g3[i,3],sens)
    
    if(lag>2){
    infetti.fino.a.qui[i]<-sum(infetti[1:10])* (pos10[i]!=0)+
                           sum(infetti[1:20])*(pos10[i]==0&&pos20[i]!=0)+
                           sum(infetti[1:30])*(pos10[i]==0&&pos20[i]==0&&pos30[i]!=0)+
                           sum(infetti[1:(lag+28)])*(pos10[i]==0&&pos20[i]==0&&pos30[i]==0)
    }else{
      infetti.fino.a.qui[i]<-sum(infetti[1:10])* (pos10[i]!=0)+
        sum(infetti[1:20])*(pos10[i]==0&&pos20[i]!=0)+
        sum(infetti[1:(lag+28)])*(pos10[i]==0&&pos20[i]==0)
    }
    
    
  }
  
  
  
  #prob_0pos<-sum(t0[1,])/Nsim
  
  
  
  #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
  prob_q10<-sum(pos10>0)/Nsim
  
  
  
  index_no10<-which(pos10==0)
  index_no20<-which(pos20==0)
  index_no30<-which(pos30==0)
  index_no10_no20<-intersect(index_no10,index_no20)
  index_no10_no20_no30<-intersect(index_no10_no20, index_no30)
  #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
  prob_q10<-sum(pos10>0)/Nsim
  #IN QUANTI CASI QUARANTENO DOPO 2 SETTIMANA
  prob_q20<-sum(pos20[index_no10]>0)/Nsim
  #IN QUANTI CASI QUARANTENO DOPO 3 SETTIMANA
  prob_q30<-sum(pos30[index_no10_no20]>0)/Nsim
  
  #IN QUANTI CASI non quaranteno 
  prob_noseen<-length(index_no10_no20_no30)/Nsim
  
  
  
  prob_q10+prob_q20+prob_q30+prob_noseen
  
  
  
  #GIORNI INFEZIONE CHE NON VEDO
  
  #  media<-sum(infetti[1:14])*prob_q14+ sum(infetti[1:7])*prob_q7+  sum(infetti[1:21])*prob_q21+sum( infetti[1:28])*prob_q28 +sum(infetti)*prob_noseen
  media=mean(infetti.fino.a.qui)
  #varianza<-(sum(infetti[1:14] )-media)^2*(prob_q14)+
  #  (sum( infetti[1:7])-media)^2*prob_q7+
  #  (sum( infetti[1:21])-media)^2*prob_q21+
  # (sum( infetti[1:28])-media)^2*prob_q28+
  # (sum(infetti)-media)^2*prob_noseen
  varianza=sd(infetti.fino.a.qui)^2
  return(list( prob_q10=prob_q10,prob_q20=prob_q20, prob_q30=prob_q30, prob_noseen=prob_noseen, media=media, varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
}


########################## 2 GRUPPI su 4 SETTIMANE G1 G2 G1 G2 

#### PIANO B1

screening_B1<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=0, lag=0 ){

# screening su met? classe a rotazione a 7,14,21,28
  
  k7.vec<-k14.vec<-k21.vec<-k28.vec<-pos7.vec<-pos14.vec<-pos21.vec<-pos28.vec<-NULL
  infetti.fino.a.qui<-NULL
  for (i in 1:Nsim){
  sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom=binom)
  infetti<-c(rep(0,lag),sir[[1]])

Inew<-c(rep(0,lag),sir$I.new) 
infetti7<-sum(Inew[1:7]) 
infetti14<-sum(Inew[1:14])
infetti21<-sum(Inew[1:21])
infetti28<-sum(Inew[1:28])
infetti_all<-sum(Inew[1:(35+lag)]) 

  
 #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato

k7<-rhyper(nn=1,m=infetti7, n=N-infetti7,k=N/2) #infetti nel primo gruppo
# infetti nel gruppo 2 alla prima settimana = infetti7-k7, quindi N/2-(infetti-k7) suscettibili

A<-rhyper(nn=1,m=infetti14-infetti7, n=N-infetti14 ,k=N/2-(infetti7-k7)) #nuovi infetti nel gruppo 2 alla settimana 2
B<-infetti14-infetti7-A #nuovi infetti nel gruppo 1 alla settimana 2
k14<-infetti7-k7+ A #infetti nel gruppo 2 alla settimana 2 

C<-rhyper(nn=1,m=infetti21-infetti14,n=N-infetti21,k=N/2-k7-B) #nuovi infetti nel gruppo 1 alla settimana 3
D<- infetti21-infetti14-C #nuovi infetti nel gruppo 2 alla settimana 3
k21<-k7+B+C #infetti nel gruppo 1 alla settimana 3

E<-rhyper(nn=1,m=infetti28-infetti21,n=N-infetti28,k=N/2-k14-D) #nuovi infetti nel gruppo 2 alla settimana 4
F<-infetti28-infetti21-E #nuovi infetti nel gruppo 1 alla settimana 4
k28<-k14+D+E #infetti nel gruppo 2 alla settimana 4


pos7<-rbinom(1,k7,sens)
pos14<-rbinom(1,k14,sens)
pos21<-rbinom(1,k21,sens)
pos28<-rbinom(1,k28,sens)

k7.vec<-c(k7.vec,k7)
k14.vec<-c(k14.vec,k14)
k21.vec<-c(k21.vec,k21)
k28.vec<-c(k28.vec,k28)
pos7.vec<-c(pos7.vec,pos7)
pos14.vec<-c(pos14.vec,pos14)
pos21.vec<-c(pos21.vec,pos21)
pos28.vec<-c(pos28.vec,pos28)

infetti.fino.a.qui[i]<-sum(infetti[1:7])*(pos7!=0)+
sum(infetti[1:14])*(pos14!=0&&pos7==0)+
sum(infetti[1:21])*(pos21!=0&&pos14==0&&pos7==0)+
sum(infetti[1:28])*(pos28!=0&&pos21==0&&pos14==0&&pos7==0)+
sum(infetti[1:(28+lag)])*(pos28==0&&pos21==0&&pos14==0&&pos7==0)

}

prob_q7<-sum(pos7.vec>0)/Nsim
prob_q14<-sum(pos7.vec==0 & pos14.vec>0)/Nsim
prob_q21<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec>0)/Nsim
prob_q28<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec==0 & pos28.vec>0)/Nsim
prob_0pos<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec==0 & pos28.vec==0)/Nsim


#media<-sum(infetti[1:14])*prob_2s+ sum(infetti[1:7])*prob_q7+ sum( infetti[1:21])*prob_0pos
media<-mean(infetti.fino.a.qui)
varianza<-sd(infetti.fino.a.qui)^2
#varianza<-(sum(infetti[1:14]) -media)^2*prob_2s+(sum(infetti[1:7])-media)^2*prob_q7+ (sum(infetti[1:21])-media)^2*prob_0pos

return(list(prob_q7=prob_q7,prob_q14=prob_q14,prob_q21=prob_q21,prob_q28=prob_q28,prob_0pos=prob_0pos, media=media, varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
}



########################## 2 GRUPPI su 4 SETTIMANE G1 G2 G1 G2 

#### PIANO B2

screening_B2<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=0, lag=0 ){
  
  # screening su met? classe a rotazione a 7,14,21,28
  
  k14.vec<-k28.vec<-pos14.vec<-pos28.vec<-NULL
  infetti.fino.a.qui<-NULL
  for (i in 1:Nsim){
    sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom=binom)
    infetti<-c(rep(0,lag),sir[[1]])
    
    Inew<-c(rep(0,lag),sir$I.new) 
       #ho messo uno perche si inizia sempre con un infetto
    infetti14<-sum(Inew[1:14])
     
    infetti28<-sum(Inew[1:28])
    infetti_all<-sum(Inew[1:(35+lag)])
    ##### QUI NON SONO SICURA CHE VADA AGGIUNTO 1
    
    
    #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato
    
    k14_g1<-rhyper(nn=1,m=infetti14, n=N-infetti14,k=N/2) #infetti nel primo gruppo a 14 gg
    k14_g2<-infetti14-k14_g1#infetti nel secondo gruppo a 14 gg  
    k14<-k14_g1
    k28_g2<-rhyper(nn=1,m=infetti28-infetti14, n=N-infetti28 ,k=N/2-k14_g2) #nuovi infetti nel gruppo 2 alla settimana 2
    k28<-k28_g2+k14_g2
      
    
   
    pos14<-rbinom(1,k14,sens)
 
    pos28<-rbinom(1,k28,sens)
    
     
    pos14.vec<-c(pos14.vec,pos14)
  
    pos28.vec<-c(pos28.vec,pos28)
    
    infetti.fino.a.qui[i]<-sum(infetti[1:14])*(pos14!=0)+
      sum(infetti[1:28])*(pos28!=0&&pos14==0)+
      sum(infetti[1:(28+lag)])*(pos28==0&&pos14==0)
  }
   
  prob_q14<-sum(  pos14.vec>0)/Nsim 
  prob_q28<-sum(  pos14.vec==0   & pos28.vec>0)/Nsim
  prob_0pos<-sum(  pos14.vec==0   & pos28.vec==0)/Nsim
  
  
  #media<-sum(infetti[1:14])*prob_2s+ sum(infetti[1:7])*prob_q7+ sum( infetti[1:21])*prob_0pos
  media<-mean(infetti.fino.a.qui)
  varianza<-sd(infetti.fino.a.qui)^2
  #varianza<-(sum(infetti[1:14]) -media)^2*prob_2s+(sum(infetti[1:7])-media)^2*prob_q7+ (sum(infetti[1:21])-media)^2*prob_0pos
  
  return(list(prob_q14=prob_q14,prob_q28=prob_q28,prob_0pos=prob_0pos, media=media, varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
}


#screening_B2()

########################## 2 GRUPPI su 4 SETTIMANE G1 G2 G1 G2 
# 
# screening_alternati_ripetuti_pooltesting<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=0, lag=0){
# 
# # screening su met? classe a rotazione a 7,14,21,28 con pooling di 6
#   
#   k7.vec<-k14.vec<-k21.vec<-k28.vec<-pos7.vec<-pos14.vec<-pos21.vec<-pos28.vec<-NULL
#   infetti.fino.a.qui<-NULL
#   for (i in 1:Nsim){
#   sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom=binom)
#   infetti<-c(rep(0,lag),sir[[1]])
#  
# Inew<-c(rep(0,lag),sir$I.new)
# infetti7<-1+sum(Inew[1:7])  #ho messo uno perche si inizia sempre con un infetto
# infetti14<-1+sum(Inew[1:14])
# infetti21<-1+sum(Inew[1:21])
# infetti28<-1+sum(Inew[1:28])
# infetti_all<-1+sum(Inew[1:(35+lag)])
# ##### QUI NON SONO SICURA CHE VADA AGGIUNTO 1
# 
#   
#  #n.b. dovremmo tener conto di quelli che guariscono ma per ora ho ignorato
# 
# k7<-rhyper(nn=1,m=infetti7, n=N-infetti7,k=N/2) #infetti nel primo gruppo
# # infetti nel gruppo 2 alla prima settimana = infetti7-k7, quindi N/2-(infetti-k7) suscettibili
# 
# A<-rhyper(nn=1,m=infetti14-infetti7, n=N-infetti14 ,k=N/2-(infetti7-k7)) #nuovi infetti nel gruppo 2 alla settimana 2
# B<-infetti14-infetti7-A #nuovi infetti nel gruppo 1 alla settimana 2
# k14<-infetti7-k7+ A #infetti nel gruppo 2 alla settimana 2 
# 
# C<-rhyper(nn=1,m=infetti21-infetti14,n=N-infetti21,k=N/2-k7-B) #nuovi infetti nel gruppo 1 alla settimana 3
# D<- infetti21-infetti14-C #nuovi infetti nel gruppo 2 alla settimana 3
# k21<-k7+B+C #infetti nel gruppo 1 alla settimana 3
# 
# E<-rhyper(nn=1,m=infetti28-infetti21,n=N-infetti28,k=N/2-k14-D) #nuovi infetti nel gruppo 2 alla settimana 4
# F<-infetti28-infetti21-E #nuovi infetti nel gruppo 1 alla settimana 4
# k28<-k14+D+E #infetti nel gruppo 2 alla settimana 4
# 
# numpos_subgroup1<-rhyper(nn=1,m=N/4,n=N/4,k7)
# numpos_subgroup2<-k7-numpos_subgroup1
# pos_subgroup1<-1*(sum(numpos_subgroup1)>0)
# pos_subgroup2<-1*(sum(numpos_subgroup2)>0)
# pos7<-rbinom(1,pos_subgroup1+pos_subgroup2,sens)
# 
# numpos_subgroup1<-rhyper(nn=1,m=N/4,n=N/4,k14)
# numpos_subgroup2<-k14-numpos_subgroup1
# pos_subgroup1<-1*(sum(numpos_subgroup1)>0)
# pos_subgroup2<-1*(sum(numpos_subgroup2)>0)
# pos14<-rbinom(1,pos_subgroup1+pos_subgroup2,sens)
# 
# numpos_subgroup1<-rhyper(nn=1,m=N/4,n=N/4,k21)
# numpos_subgroup2<-k21-numpos_subgroup1
# pos_subgroup1<-1*(sum(numpos_subgroup1)>0)
# pos_subgroup2<-1*(sum(numpos_subgroup2)>0)
# pos21<-rbinom(1,pos_subgroup1+pos_subgroup2,sens)
# 
# numpos_subgroup1<-rhyper(nn=1,m=N/4,n=N/4,k28)
# numpos_subgroup2<-k28-numpos_subgroup1
# pos_subgroup1<-1*(sum(numpos_subgroup1)>0)
# pos_subgroup2<-1*(sum(numpos_subgroup2)>0)
# pos28<-rbinom(1,pos_subgroup1+pos_subgroup2,sens)
# 
# k7.vec<-c(k7.vec,k7)
# k14.vec<-c(k14.vec,k14)
# k21.vec<-c(k21.vec,k21)
# k28.vec<-c(k28.vec,k28)
# pos7.vec<-c(pos7.vec,pos7)
# pos14.vec<-c(pos14.vec,pos14)
# pos21.vec<-c(pos21.vec,pos21)
# pos28.vec<-c(pos28.vec,pos28)
# 
# infetti.fino.a.qui[i]<-sum(infetti[1:7])*(pos7!=0)+
# sum(infetti[1:14])*(pos14!=0&&pos7==0)+
# sum(infetti[1:21])*(pos21!=0&&pos14==0&&pos7==0)+
# sum(infetti[1:28])*(pos28!=0&&pos21==0&&pos14==0&&pos7==0)+
# sum(infetti[1:(35+lag)])*(pos28==0&&pos21==0&&pos14==0&&pos7==0)
# 
# }
# 
# prob_q7<-sum(pos7.vec>0)/Nsim
# prob_q14<-sum(pos7.vec==0 & pos14.vec>0)/Nsim
# prob_q21<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec>0)/Nsim
# prob_q28<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec==0 & pos28.vec>0)/Nsim
# prob_0pos<-sum(pos7.vec==0 & pos14.vec==0 & pos21.vec==0 & pos28.vec==0)/Nsim
# 
# 
# #media<-sum(infetti[1:14])*prob_2s+ sum(infetti[1:7])*prob_q7+ sum( infetti[1:21])*prob_0pos
# media<-mean(infetti.fino.a.qui)
# varianza<-sd(infetti.fino.a.qui)^2
# #varianza<-(sum(infetti[1:14]) -media)^2*prob_2s+(sum(infetti[1:7])-media)^2*prob_q7+ (sum(infetti[1:21])-media)^2*prob_0pos
# 
# return(list(prob_q7=prob_q7,prob_q14=prob_q14,prob_q21=prob_q21,prob_q28=prob_q28,prob_0pos=prob_0pos, media=media, varianza=varianza))
# }
# 


#########################

##### PIANO A2

screening_A2<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=1, lag=0){

# screening su tutti a 14 e 28

  infetti.fino.a.qui<-NULL
  
  pos28T.vec<-NULL
  pos14T.vec<-NULL
  for(j in 1:Nsim){
    sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom)
 
    infetti<-c(rep(0,lag),sir[[1]])
    Inew<-c(rep(0,lag),sir$I.new)
    
    infetti14<-sum(Inew[1:14])
    infetti_all<-sum(Inew[1:(35+lag)])
    pos14<-rbinom(1,infetti14,sens)
    pos14T.vec<-c(pos14T.vec,pos14)
    infetti28<-sum(Inew[1:28])
    pos28<-rbinom(1,infetti28,sens)
    pos28T.vec<-c(pos28T.vec,pos28)
    
    infetti.fino.a.qui[j]<-sum(infetti[1:14])*(pos14!=0)+sum(infetti[1:28])*(pos14==0&&pos28!=0)+sum(infetti[1:(28+lag)])*(pos14==0&&pos28==0)
  }
  
   
  index_no14<-which(pos14T.vec==0)
  index_no28<-which(pos28T.vec==0)
  index_no28_no14<-intersect(index_no14,index_no28)  
  #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
    #IN QUANTI CASI QUARANTENO DOPO 2 SETTIMANA
  prob_q14<-sum(pos14T.vec>0)/Nsim
  prob_q28<-sum(pos28T.vec[index_no14]>0)/Nsim
  prob_0pos<-length(index_no28_no14)/Nsim
  
  media<-mean(infetti.fino.a.qui)
varianza<-sd(infetti.fino.a.qui)^2  
  #varianza<-(sum(round(infetti[1:14],0))-media)^2*(1-prob_0pos)+ (sum(round(infetti[1:21],0))-media)^2*(1-prob_0pos)
  return(list(pos14T.vec=pos14T.vec,prob_q14=prob_q14, prob_q28=prob_q28,  prob_0pos=prob_0pos,media=media , varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
}
 
#########################

##### PIANO A1

screening_A1<- function(N=20, sens=0.7, R0=3, tempo.sim=60, T.exit=14, Nsim=10000, binom=1, lag=0){
 

  infetti.fino.a.qui<-NULL

  pos7T.vec<-pos21T.vec<-pos28T.vec<-pos14T.vec<-NULL
  for(j in 1:Nsim){
    sir<-SIR(1,N,0,R0,T.exit,tempo.sim, binom)
    infetti<-c(rep(0,lag),sir[[1]])
    Inew<-c(rep(0,lag),sir$I.new)
    infetti7<-sum(Inew[1:7])
    pos7<-rbinom(1,infetti7,sens)
    pos7T.vec<-c(pos7T.vec,pos7)
    infetti14<-sum(Inew[1:14])
    pos14<-rbinom(1,infetti14,sens)
    pos14T.vec<-c(pos14T.vec,pos14)
    infetti21<-sum(Inew[1:21])
    pos21<-rbinom(1,infetti21,sens)
    pos21T.vec<-c(pos21T.vec,pos21)
    infetti28<-sum(Inew[1:28])
    pos28<-rbinom(1,infetti28,sens)
    pos28T.vec<-c(pos28T.vec,pos28)
    infetti_all<-sum(Inew[1:(35+lag)])
    
    infetti.fino.a.qui[j]<-sum(infetti[1:7])*(pos7!=0)+sum(infetti[1:14])*(pos14!=0&&pos7==0)+sum(infetti[1:21])*(pos21!=0&&pos7==0&&pos14==0)+
    sum(infetti[1:28])*(pos7==0&&pos14==0&&pos21==0&&pos28!=0)+ sum(infetti[1:(28+lag)])*(pos7==0&&pos14==0&&pos21==0&&pos28==0)
  }
  
   
  index_no7<-which(pos7T.vec==0)
  index_no14<-which(pos14T.vec==0)
  index_no21<-which(pos21T.vec==0)
  index_no28<-which(pos28T.vec==0)
  index_no7_no14<-intersect(index_no7,index_no14)
  index_no7_no14_no21<-intersect(index_no7_no14, index_no21)  
  index_no7_no14_no21_no28<-intersect(index_no7_no14_no21, index_no28) 
  #IN QUANTI CASI QUARANTENO DOPO 1 SETTIMANA
    #IN QUANTI CASI QUARANTENO DOPO 2 SETTIMANA
  prob_q7<-sum(pos7T.vec>0)/Nsim
  prob_q14<-sum(pos14T.vec[index_no7]>0)/Nsim
  prob_q21<-sum(pos21T.vec[index_no7_no14]>0)/Nsim
  prob_q28<-sum(pos28T.vec[index_no7_no14_no21]>0)/Nsim
  prob_0pos<-length(index_no7_no14_no21_no28)/Nsim

  
  media<-mean(infetti.fino.a.qui)
varianza<-sd(infetti.fino.a.qui)^2  
  #varianza<-(sum(round(infetti[1:14],0))-media)^2*(1-prob_0pos)+ (sum(round(infetti[1:21],0))-media)^2*(1-prob_0pos)
  return(list(pos14T.vec=pos14T.vec,prob_q7=prob_q7,prob_q14=prob_q14, prob_q21=prob_q21, prob_q28=prob_q28,  prob_0pos=prob_0pos,media=media , varianza=varianza, infetti.fino.a.qui=infetti.fino.a.qui))
}
 