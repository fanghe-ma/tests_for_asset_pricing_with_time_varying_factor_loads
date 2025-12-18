* Galvao, Montes-Rojas, Olmo
* code for tests
* It requires *.csv matrices
* It produces the test results

#R-function adapted from Ando and Bai (2015)
rm(list=ls())

#T # length of time series
#N # number of indiviuals
#k # dimension of explanatory variables (including constant)
#r #dimension of common factors
#AY # Panel of dependent variable, one column for each individual
#AX # Explanatory variables, k columns for each individual

AndoBai<- function(AY,AX,AXb,Av,r){

N<-length(AY[1,])
T<-length(AY[,1])
k<-length(AX[1,])/N

Y <- AY

PredXB <- matrix(0,nrow=T,ncol=sum(N))
Bes <- 0*matrix(runif(k*sum(N),-1,1),ncol=sum(N))

for(j in 1:sum(N)){
X <- AX[,(k*(j-1)+1):(k*j)]; y <- Y[,j]
fit <- lm(y~X-1)
Bes[,j] <- (fit$coefficients)
PredXB[,j] <- X%*%as.vector(Bes[,j])
}

Bes_noF<-Bes
u<-AY-PredXB

Y <- AY-PredXB

VEC <- eigen(Y%*%t(Y))$vectors
F <- sqrt(T)*(VEC)[,1:r]
L <- t(t(F)%*%Y/T)

FG <- F; LG <- L

PredG <- matrix(0,nrow=T,ncol=sum(N))
for(j in 1:sum(N)){PredG[,j] <- FG%*%LG[j,]}


for(ite in 1:100){
Bes.old <- Bes
Y <- AY-PredG
for(j in 1:sum(N)){
X <- AX[,(k*(j-1)+1):(k*j)]; y <- Y[,j]
fit <- lm(y~X-1)
Bes[,j] <- (fit$coefficients)
PredXB[,j] <- X%*%as.vector(Bes[,j])
}
Y <- AY-PredXB
VEC <- eigen(Y%*%t(Y))$vectors
F <- sqrt(T)*(VEC)[,1:r]
L <- t(t(F)%*%Y/T)
FG <- F; LG <- L
PredG <- matrix(0,nrow=T,ncol=sum(N))
for(j in 1:sum(N)){PredG[,j] <- FG%*%LG[j,]}
if(sum( abs(Bes.old-Bes) )<=10^-9){break}
}


#---------------Swamy test

H0Bes <- Bes%*%rep(1,len=N)/N
H0Bes[1,1] <-0 

H0Bes_noF <- Bes_noF%*%rep(1,len=N)/N
H0Bes_noF[1,1] <-0 

ER <- AY-PredXB-PredG

XX <- S <- L <- Omega_gr <- Omega_ngr <- Omega_gr_noF <- Omega_ngr_noF <- matrix(0,k*N,k*N)
#Eta <- matrix(0,k*N,1)

M <- diag(1,T)-FG%*%solve(t(FG)%*%FG)%*%t(FG)

BBB <- BBB_noF <- matrix(0,k*N,1)

for(j in 1:N){BBB[((j-1)*k+1):(j*k),1] <- Bes[,j]-H0Bes}
for(j in 1:N){BBB_noF[((j-1)*k+1):(j*k),1] <- Bes_noF[,j]-H0Bes_noF}

BBB<-BBB
BBB_noF<-BBB_noF
Bes<-Bes
Bes_noF<-Bes_noF

s <- sum( ER^2 )/(T*N-N*k-(T+N)*r)

for(j in 1:N){

X <- AX[,(k*(j-1)+1):(k*j)]; 
Xb <- as.matrix(AXb[,((k-1)*(j-1)+1):((k-1)*j)],nrow(AXb),k-1); 
v <- as.matrix(Av[,((k-1)*(j-1)+1):((k-1)*j)],nrow(AXb),k-1);

#From Stata
#mat factoromega=invsym(XX)*XvXv*invsym(XX)
#* To form term1=Xhat_eta'X*factoromega*X'Xhat_eta
#	mat term1=invsym(XhatXhat)*(Xhat_eta'*X*factoromega*X'*Xhat_eta)*invsym(XhatXhat)

aux1<-s*t(X)%*%M%*%X/T
aux2<-aux1*0
MX<-M%*%X
for(h in 2:k){
	aux2[h,h]<-1/T*Bes[h,j]^2*(sum(MX[,h]*Xb[,h-1]))^2*sum(Xb[,h-1]^2*v[,h-1]^2)/((sum(Xb[,h-1]^2))^2)
}

Omega_gr[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <- aux1+aux2
Omega_ngr[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <- aux1

aux1<-sum(u[,j]^2)/(T-k)*t(X)%*%X/T

Omega_gr_noF[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <- aux1+aux2
Omega_ngr_noF[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <- aux1

XX[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <-  t(X)%*%X/T

S[((j-1)*k+1):(j*k),((j-1)*k+1):(j*k)] <- t(X)%*%M%*%X/T
#Eta[((j-1)*k+1):(j*k),1] <- t(X)%*%M%*%ER[,j]/sqrt(T)

for(l in 1:N){
Xl <- AX[,(k*(l-1)+1):(k*l)]; 
a <- sum( t(LG[j,])%*%solve(t(LG)%*%LG/N)%*%LG[l,] )
L[((j-1)*k+1):(j*k),((l-1)*k+1):(l*k)] <- a*t(X)%*%M%*%Xl/T
}

}

V_ab_gr <<- (S-L/N)%*%solve(Omega_gr)%*%t(S-L/N)
V_ab_ngr <<- (S-L/N)%*%solve(Omega_ngr)%*%t(S-L/N)
V_ab_gr_noF <- XX%*%solve(Omega_gr_noF)%*%t(XX)
V_ab_ngr_noF <- XX%*%solve(Omega_ngr_noF)%*%t(XX)

Swamy_ab_gr <- T*t(BBB)%*%V_ab_gr%*%BBB
Swamy_ab_ngr <- T*t(BBB)%*%V_ab_ngr%*%BBB

Swamy_ab_gr_noF <- T*t(BBB_noF)%*%V_ab_gr_noF%*%BBB_noF
Swamy_ab_ngr_noF <- T*t(BBB_noF)%*%V_ab_ngr_noF%*%BBB_noF

vector_a<-BBB*0
for(j in 1:N){
	vector_a[(j-1)*k+1,1]<-1
}
vector_b<-BBB^0-vector_a


Swamy_a_gr <- T*t(BBB*vector_a)%*%V_ab_gr%*%(BBB*vector_a)
Swamy_a_ngr <- T*t(BBB*vector_a)%*%V_ab_ngr%*%(BBB*vector_a)
Swamy_a_gr_noF <- T*t(BBB_noF*vector_a)%*%V_ab_gr_noF%*%(BBB_noF*vector_a)
Swamy_a_ngr_noF <- T*t(BBB_noF*vector_a)%*%V_ab_ngr_noF%*%(BBB_noF*vector_a)

Swamy_b_gr <- T*t(BBB*vector_b)%*%V_ab_gr%*%(BBB*vector_b)
Swamy_b_ngr <- T*t(BBB*vector_b)%*%V_ab_ngr%*%(BBB*vector_b)
Swamy_b_gr_noF <- T*t(BBB_noF*vector_b)%*%V_ab_gr_noF%*%(BBB_noF*vector_b)
Swamy_b_ngr_noF <- T*t(BBB_noF*vector_b)%*%V_ab_ngr_noF%*%(BBB_noF*vector_b)

Swamy_ab_gr <- as.numeric( (Swamy_ab_gr-(N-1)*(k-1)-N)/sqrt(2*((N-1)*(k-1)+N)) )
Swamy_ab_ngr <- as.numeric( (Swamy_ab_ngr-(N-1)*(k-1)-N)/sqrt(2*((N-1)*(k-1)+N)) )
Swamy_ab_gr_noF <- as.numeric( (Swamy_ab_gr_noF-(N-1)*(k-1)-N)/sqrt(2*((N-1)*(k-1)+N)) )
Swamy_ab_ngr_noF <- as.numeric( (Swamy_ab_ngr_noF-(N-1)*(k-1)-N)/sqrt(2*((N-1)*(k-1)+N)) )

Swamy_a_gr <- as.numeric( (Swamy_a_gr-N)/sqrt(2*N) )
Swamy_a_ngr <- as.numeric( (Swamy_a_ngr-N)/sqrt(2*N) )
Swamy_a_gr_noF <- as.numeric( (Swamy_a_gr_noF-N)/sqrt(2*N) )
Swamy_a_ngr_noF <- as.numeric( (Swamy_a_ngr_noF-N)/sqrt(2*N) )

Swamy_b_gr <- as.numeric( (Swamy_b_gr-(N-1)*(k-1))/sqrt(2*((N-1)*(k-1))) )
Swamy_b_ngr <- as.numeric( (Swamy_b_ngr-(N-1)*(k-1))/sqrt(2*((N-1)*(k-1))) )
Swamy_b_gr_noF <- as.numeric( (Swamy_b_gr_noF-(N-1)*(k-1))/sqrt(2*((N-1)*(k-1))) )
Swamy_b_ngr_noF <- as.numeric( (Swamy_b_ngr_noF-(N-1)*(k-1))/sqrt(2*((N-1)*(k-1))) )

c(N,T,r,k,
Swamy_a_gr,Swamy_a_ngr,
Swamy_a_gr_noF,Swamy_a_ngr_noF,
Swamy_b_gr,Swamy_b_ngr,
Swamy_b_gr_noF,Swamy_b_ngr_noF,
Swamy_ab_gr,Swamy_ab_ngr,
Swamy_ab_gr_noF,Swamy_ab_ngr_noF)

}

setwd("F:/SDF/realized data full sample")


AX_1_2014<-as.matrix(read.table("AX_1_2014.csv", header=FALSE,     sep=","))
AXb_1_2014<-as.matrix(read.table("AXb_1_2014.csv", header=FALSE,     sep=","))
Av_1_2014<-as.matrix(read.table("Av_1_2014.csv", header=FALSE,     sep=","))
AY1_1_2014<-as.matrix(read.table("AY1_1_2014.csv", header=FALSE,     sep=","))
AY2_1_2014<-as.matrix(read.table("AY2_1_2014.csv", header=FALSE,     sep=","))
AX_3_2014<-as.matrix(read.table("AX_3_2014.csv", header=FALSE,     sep=","))
AXb_3_2014<-as.matrix(read.table("AXb_3_2014.csv", header=FALSE,     sep=","))
Av_3_2014<-as.matrix(read.table("Av_3_2014.csv", header=FALSE,     sep=","))
AY1_3_2014<-as.matrix(read.table("AY1_3_2014.csv", header=FALSE,     sep=","))
AY2_3_2014<-as.matrix(read.table("AY2_3_2014.csv", header=FALSE,     sep=","))
AX_5_2014<-as.matrix(read.table("AX_5_2014.csv", header=FALSE,     sep=","))
AXb_5_2014<-as.matrix(read.table("AXb_5_2014.csv", header=FALSE,     sep=","))
Av_5_2014<-as.matrix(read.table("Av_5_2014.csv", header=FALSE,     sep=","))
AY1_5_2014<-as.matrix(read.table("AY1_5_2014.csv", header=FALSE,     sep=","))
AY2_5_2014<-as.matrix(read.table("AY2_5_2014.csv", header=FALSE,     sep=","))

r=2
ini=1
end=204
AndoBai(AY1_1_2014[ini:end,],AX_1_2014[ini:end,],AXb_1_2014[ini:end,],Av_1_2014[ini:end,],r)
AndoBai(AY1_3_2014[ini:end,],AX_3_2014[ini:end,],AXb_3_2014[ini:end,],Av_3_2014[ini:end,],r)
AndoBai(AY1_5_2014[ini:end,],AX_5_2014[ini:end,],AXb_5_2014[ini:end,],Av_5_2014[ini:end,],r)
ini=1
end=120
AndoBai(AY1_1_2014[ini:end,],AX_1_2014[ini:end,],AXb_1_2014[ini:end,],Av_1_2014[ini:end,],r)
AndoBai(AY1_3_2014[ini:end,],AX_3_2014[ini:end,],AXb_3_2014[ini:end,],Av_3_2014[ini:end,],r)
AndoBai(AY1_5_2014[ini:end,],AX_5_2014[ini:end,],AXb_5_2014[ini:end,],Av_5_2014[ini:end,],r)
ini=40
end=160
AndoBai(AY1_1_2014[ini:end,],AX_1_2014[ini:end,],AXb_1_2014[ini:end,],Av_1_2014[ini:end,],r)
AndoBai(AY1_3_2014[ini:end,],AX_3_2014[ini:end,],AXb_3_2014[ini:end,],Av_3_2014[ini:end,],r)
AndoBai(AY1_5_2014[ini:end,],AX_5_2014[ini:end,],AXb_5_2014[ini:end,],Av_5_2014[ini:end,],r)
ini=80
end=204
AndoBai(AY1_1_2014[ini:end,],AX_1_2014[ini:end,],AXb_1_2014[ini:end,],Av_1_2014[ini:end,],r)
AndoBai(AY1_3_2014[ini:end,],AX_3_2014[ini:end,],AXb_3_2014[ini:end,],Av_3_2014[ini:end,],r)
AndoBai(AY1_5_2014[ini:end,],AX_5_2014[ini:end,],AXb_5_2014[ini:end,],Av_5_2014[ini:end,],r)

