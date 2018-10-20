#' A function to calculate the deterministic estimates of the number
#' infected versus time and the average number of descendant infections (ANDI)
#' for a deterministic Susceptible, Infected, Recovered (SIR) model
#' @export
SIR_solve_model = function(N,
                           R0,
                           fsusc,
                           delta_t=0.01){

  I_0 = 1      
  S_0 = fsusc*N-I_0
  R_0 = (1-fsusc)*N
 
  tend = 1*365

  # note that time is in units of 1/gamma
  tbeg  = 0           
  gamma = 1/1         
  beta = R0*gamma     

  vparameters = c(gamma=gamma,beta=beta)
  inits       = c(S=S_0,I=I_0,R=R_0)

  iend = 1
  while (iend>0.0001){
    vt = seq(0,tend,delta_t)
    sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))
    sirmodel = subset(sirmodel,!is.na(S+I+R))
    iend = sirmodel$I[nrow(sirmodel)]
    tend = tend*2
  }
  minS = min(sirmodel$S,na.rm=T)
  final = (max(sirmodel$S,na.rm=T)-min(sirmodel$S,na.rm=T))/N

  ################################################################################## 
  ################################################################################## 
  # there is a problem when I goes to less than 1, but S-Sinfty is still non-negligble
  # For high R0, (S-Sinfty) drops quickly and is close to zero when I drops below 1.
  # In contrast, for R0=1.05 (S-Sinfty) is still over 16 when I
  # drops below 1
  ################################################################################## 
  a = sirmodel
  a$newI = beta*a$S*a$I/N
  a$wprob = a$newI/sum(a$newI)
  a$wfraction_responsible_for = 1/(a$I)
  a$wfraction_responsible_for[a$wfraction_responsible_for>1] = 0
  iind = which(a$wfraction_responsible_for==0)

  n = iind[1]-1
  a$wprob[n] = sum(a$wprob[iind])
  a$Nafter = a$S-minS
  a$NDI = a$Nafter*a$wfraction_responsible_for

  ANDI = weighted.mean(a$NDI,a$wprob)

  return(list(results    = a,
              final_size = final,
              ANDI       = ANDI
             ))
} # end solve_model function definition



