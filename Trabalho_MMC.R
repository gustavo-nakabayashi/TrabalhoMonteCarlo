library(triangle)
M <- 10000 #Numero de iterações
Vs <- c(2.0361,6.0881,10.0229,14.048,19.035) #Medidas do Padrao
faixa_Vs <- 100 #Escala de 100v do Multimetro padrao
delta_Vs <- matrix(c(1:5*M), nrow=M, ncol=5, byrow = TRUE) #Incerteza do Padrao
delta_Vs_sd <-c(1:5) #Desvio padrão da medida do Padrao
Vi <- matrix(c(1:5*M), nrow=M, ncol=5, byrow = TRUE) #Distribuição dos valores medidos
delta_Vi <- matrix(c(1:5*M), nrow=M, ncol=5, byrow = TRUE) #Erro da resolução do multimetro
delta_Vi_sd <- c(1:5)#Desvio padrao do erro da resolução do multimetro
Ex <- matrix(c(1:5*M), nrow=M, ncol=5, byrow = TRUE) #Erro da medida
Ex_media <- c(1:5)#Media do erro
Ex_sd <- c(1:5)#Media do erro
Medidas = matrix(c(2.02, 2.01, 2.01,2.00,2.01,
                   6.07, 6.07,6.07, 6.07,6.06,
                   10.00, 10.00, 9.99, 9.99, 10.00,
                   14.02, 14.01, 14.00, 14.00, 14.01,
                   18.99, 18.99, 19.00, 18.99, 18.99), 
                 nrow=5, ncol=5, byrow = TRUE) 
Medidas_sd <-c(1:5)#Desvio padrao das medidas
Medidas_mean<-c(1:5)#Valor medio das medidas
LPU_distr <- matrix(c(1:5*M), nrow=M, ncol=5, byrow = TRUE)
Incerteza_LPU_095 <- -c(1:5)#Valor medio das medidas
lim_Sup_LPU <- c(1:5)
lim_Inf_LPU <- c(1:5)
Uc <- c(1:5)
Veff <-c(1:5)
for(i in 1:5)
{
  delta_Vs_sd[i] <- 0.0045*Vs[i]/100 + 0.0006*faixa_Vs/100 #±( % of reading + % of range )
  delta_Vs[,i] <- rnorm(M, 0, delta_Vs_sd[i])
  delta_Vi_sd[i] <- (0.01)/sqrt(12)
  Medidas_sd[i] <- sd(Medidas[i,])
  Medidas_mean[i] <- mean(Medidas[i,])
  delta_Vi[,i] <- runif(n=M, min=-0.005, max=0.005)
  Vi[,i] <- propagate:::rst(M, mean = Medidas_mean[i], sd = Medidas_sd[i], df = 4)
  #print(delta_Vs[i])
}


for(i in 1:5)
{
  for(j in 1:M)
  {
    Ex[j,i] <- Vi[j,i] + delta_Vi[j,i] - Vs[i] - delta_Vs[j,i]
  }
}

print("Os erros das medidas Ex são respectivamente")
for(i in 1:5)
{
  Ex_media[i] <- mean(Ex[,i])
  Ex_sd[i] <- sd(Ex[,i])
  print(format(round(Ex_media[i], 2), nsmall = 2))
}

######## LPU
for(i in 1:5)
{
  Uc[i] <- sqrt(Medidas_sd[i]^2 + delta_Vs_sd[i]^2 + delta_Vi_sd[i]^2)
  LPU_distr[,i] <- propagate:::rst(M, mean = (Medidas_mean[i] - Vs[i]), sd = Uc[i], df = 4)
  lim_Sup_LPU[i] <- quantile(LPU_distr[,i], 0.975)
  lim_Inf_LPU[i] <- quantile(LPU_distr[,i], 0.025)
  Incerteza_LPU_095[i] <- (lim_Sup_LPU[i]-lim_Inf_LPU[i])/2
  print(paste0("Incerteza Expandida LPU:",i," ", format(round(Incerteza_LPU_095[i], 2), nsmall = 2)))
  Veff[i] <- (Uc[i]^4)/((Medidas_sd[i]^4)/4)
  print(paste0("Veff LPU:",i," ", format(round(Veff[i], 2), nsmall = 2)))
  
}
#########


for (i in 1:5) 
{
  Ex[,i] <- sort(Ex[,i])
}


lim_Sup <- c(1:5)
lim_Inf <- c(1:5)
Veff_MMC <- c(1:5)
r<-M*(1-0.95)/2
q<- 0.95*M 
for(i in 1:5)
{
  
  lim_Sup[i] <- (Ex[(r+q),i])
  lim_Inf[i] <- (Ex[(r),i])
  # print(paste0("Limite superior é igual a: ",lim_Sup[i]))
  # print(paste0("Limite inferior é igual a: ",lim_Inf[i]))
  Incerteza[i] <- ( lim_Sup[i] - lim_Inf[i])/2
  # print(paste0("Incerteza MMC:",i," ", format(round(Incerteza[i], 2), nsmall = 2)))
  print(paste0("Limite inferior é: ", format(round(lim_Inf[i], 5), nsmall = 2), "Limite Superior é: ", format(round(lim_Sup[i], 5), nsmall = 2)))
  print(paste0("Limite inferior é: ", format(round(lim_Inf_LPU[i], 5), nsmall = 2), "Limite Superior é: ", format(round(lim_Sup_LPU[i], 5), nsmall = 2)))
  
  plot(ecdf(Ex[,i]), xlim=c(-.08,.01),col = "red",lwd=2, main=paste0("CDF de Ex ",round(Vs[i]),"V"), ylab="Probabilidade (Ex)", xlab="Ex(V)")
  abline(v = mean(Ex[,i]), untf = FALSE)
  abline(v = lim_Sup[i], untf = FALSE,col = "red",lwd=2)
  abline(v = lim_Inf[i], untf = FALSE,col = "red",lwd=2)
  
  plot(density(Ex[,i]), xlim=c(-.08,.01),ylim=c(.0,90),col = "red",lwd=2, main=paste0("CDF de Ex ",round(Vs[i]),"V"),
       ylab="Probabilidade (Ex)",
       xlab="Ex(V)")
  lines(density(LPU_distr[,i]),col = "blue",lwd=2)
  legend(-0.013, 55, legend=c("MMC", "LPU"),
        col=c("red", "blue"), lty=1:1, cex=0.8,
        title="Métodos", text.font=4, bg='lightblue')
  abline(v = mean(Ex[,i]), untf = FALSE)
  abline(v = lim_Sup[i], untf = FALSE,col = "red",lwd=2)
  abline(v = lim_Inf[i], untf = FALSE,col = "red",lwd=2)
  abline(v = lim_Sup_LPU[i], untf = FALSE,col = "blue",lwd=2)
  abline(v = lim_Inf_LPU[i], untf = FALSE,col = "blue",lwd=2)
  Veff_MMC[i] <- (Ex_sd[i]^4)/((Medidas_sd[i]^4)/4)
  print(paste0("Veff MMC:",i," ", format(round(Veff_MMC[i], 2), nsmall = 2)))
}

int_A <- c(M*0.95):M#Interlavo de abrangência
qE <- c(M*0.95):M#Quantil esquerdo
k<- 1 
menor <- 100000
for (i in (M*0.95):M)
{
  int_A[k] = ((k-1)*1/M)
  qE[k] = Ex[i,2] - Ex[k,2]
  
  if (qE[k] < menor)
  {
    menor <- qE[k]
  }
  k <- k+1
}
plot(int_A,qE,col = "red",lwd=2, type = "l", ylab="Intervalo de abrangencia",
     xlab="Quantil Esquerdo", main="Intervalo de abrangência x quantil esquerdo 6V")

