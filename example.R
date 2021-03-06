
rm(list = ls(all = TRUE))

#require("devtools")
#install_github("smtowers/ANDI")
require("ANDI")
require("deSolve")

N = 1000
I0 = 1
S0 = N-I0
gamma = 1
R0 = 1.5

set.seed(78998026)
a = SIR_agent(N,I0,S0,gamma,R0)
b = SIR_solve_model(N,R0,S0/N,0.1)

par(mfrow=c(1,1))
plot(a$time,a$I,xlab="Time, in units of 1/gamma",ylab="Prevalence",main=paste("R0=",R0," and N=",N,sep=""))
lines(b$results$time,b$results$I,col=2,lwd=4)
legend("topright",legend=c("ABMC","deterministic"),col=c(1,2),bty="n",lwd=4,cex=0.8)

cat("Estimated ANDI, ABMC and deterministic: ",round(mean(a$ANDI),1),"  ",round(b$ANDI,1),"\n",sep="")
cat("Estimated final size, ABMC and deterministic: ",a$final_size,"  ",b$final_size,"\n",sep="")

