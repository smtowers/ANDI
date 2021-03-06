% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIR_agent.R
\name{SIR_agent}
\alias{SIR_agent}
\title{An R function to perform an Agent Based Monte Carlo simulation of a
Susceptible, Infected, Recovered (SIR) disease model with homogeneous mixing
of the population.  The function keeps track of who infected whom, and
calculates the average number of descendant infections (ANDI) for the 
individuals infected in the outbreak}
\usage{
SIR_agent(N, I_0, S_0, gamma, R0, delta_t = 0)
}
\arguments{
\item{N}{The population size}

\item{I_0}{The initial number infected}

\item{S_0}{The initial number suscetible}

\item{gamma}{The recovery rate}

\item{R0}{The reproduction number}

\item{delta_t}{time step to solve model (if 0 then dynamically calculated time step used)}
}
\value{
list list of model compartment values vs time, and ANDI
}
\description{
An R function to perform an Agent Based Monte Carlo simulation of a
Susceptible, Infected, Recovered (SIR) disease model with homogeneous mixing
of the population.  The function keeps track of who infected whom, and
calculates the average number of descendant infections (ANDI) for the 
individuals infected in the outbreak
}
\examples{
\dontrun{
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
  plot(a$time,a$I,xlab="Time",ylab="Prevalence",main=paste("R0=",R0," and N=",N,sep=""))
  lines(b$results$time,b$results$I,col=2,lwd=4)
  legend("topright",legend=c("ABMC","deterministic"),col=c(1,2),bty="n",lwd=4,cex=0.8)
  
  cat("Est ANDI, ABMC and deterministic: ",round(mean(a$ANDI),1),"  ",round(b$ANDI,1),"\\n",sep="")
  cat("Est final size, ABMC and deterministic: ",a$final_size,"  ",b$final_size,"\\n",sep="")
}

}
