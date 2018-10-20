#' These are the differential equations for a Susceptibe, Infected, Recovered (SIR)
#' determinstic model, input to the lsoda() method in the deSolve library
#' @export
SIRfunc=function(t
                ,x
                ,vparameters
                ){
   S = x[1]  
   I = x[2]  
   R = x[3]  

   S[S<0] = 0
   I[I<0] = 0
   R[R<0] = 0
   I[is.na(I)] = 0

   with(as.list(vparameters),{
      npop=S+I+R
      dS = -beta*S*I/npop
      dI = +beta*S*I/npop - gamma*I
      dR = +gamma*I
      out = c(dS,dI,dR)
      list(out)
   })
} # end SIRfunc function definition


