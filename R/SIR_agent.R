#' An R function to perform an Agent Based Monte Carlo simulation of a
#' Susceptible, Infected, Recovered (SIR) disease model with homogeneous mixing
#' of the population.  The function keeps track of who infected whom, and
#' calculates the average number of descendant infections (ANDI) for the 
#' individuals infected in the outbreak
#' @export
SIR_agent = function(N         # population size
                    ,I_0       # initial number infected
                    ,S_0       # initial number susceptible
                    ,gamma     # recovery rate in days^{-1}
                    ,R0        # reproduction number
                    ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
                    ){

   ########################################################################
   # begin by calculating the transmission rate, beta, from
   # R0 and gamma
   #
   # beta is the number of contacts sufficient to transmit infection per unit time
   #
   # Thus, in time delta_t, a susceptible person will contact a Poisson random
   # number of infected people, with mean beta*I*delta_t/N
   #
   # The probability an infected person will recover in time
   # delta_t is p=1-exp(-gamma*delta_t)
   # (note that p does not depend on t!  This is due to the 
   #  memoryless nature of the Exponential distribution)
   ########################################################################
   beta = R0*gamma 

   ########################################################################
   # now set up the state vector of the population
   # vstate = 0   means susceptible
   # vstate = 1   means infected    
   # vstate = 2   means recovered   
   # randomly pick I_0 of the people to be infected
   ########################################################################
   vstate = rep(2,N)
   vstate[1:S_0] = 0
   index_inf = (N-I_0+1):N
   vstate[index_inf] = 1     # these are the infected people
   zwho  = rep(0,N)
   ztime = rep(-1,N)
   z_linfected = rep(0,N)
   znum_descendant_infections = rep(0,N)
   zlist = list()
   zgeneration = list()
   for (i in 1:N) zlist[[i]] = numeric(0)
   for (i in 1:N) zgeneration[[i]] = numeric(0)

   ########################################################################
   # now begin the time steps
   ########################################################################
   t     = 0
   S     = S_0
   I     = I_0
   vS    = numeric(0)  # this will be filled with the number of susceptibles in the population over time
   vI    = numeric(0)  # this will be filled with the number of infectives in the population over time
   vtime = numeric(0)  # this will be filled with the time steps

   while (length(vstate[vstate==1])>0){ # continue the simulation until we have no more infectious people
      S = length(vstate[vstate==0])  # count the number of susceptibles, based on the state vector 
      I = length(vstate[vstate==1])  # count the number of infectives, based on the state vector

      vS = append(vS,S)  # append this info to vectors that we will return from the function
      vI = append(vI,I)

      vtime = append(vtime,t)

      #cat(t,S,I,"\n")

      ########################################################################
      # calculate the dynamical time step for the simulation.  Otherwise
      # use a fixed time step
      ########################################################################
      deltat=delta_t
      if (delta_t==0) deltat = 1/(beta*I*S/N + gamma*I)
      recover_prob = (1-exp(-gamma*deltat)) # the probability an infective recovers in this time step

      ########################################################################
      # sample Poisson distributed random numbers to simulated the
      # number of infected people contacted by each person
      ########################################################################
      avg_num_infected_people_contacted = beta*I*deltat/N
      vnum_infected_people_contacted = rpois(N,avg_num_infected_people_contacted) 
      vprob = runif(N)   # sample uniform random numbers
      vnewstate = vstate # copy the state vector to a temporary vector used for calculations

      ########################################################################
      # Infected people recover if the sampled uniform random
      # number is less than the recovery probability
      ########################################################################
      vnewstate[vstate==1&vprob<recover_prob] = 2              

      ########################################################################
      # if the person was infectious at any point in time, mark them
      # as one of the infected people in the outbreak
      ########################################################################
      l = which(vstate==1)
      z_linfected[l] = 1

      ########################################################################
      # If a susceptible contacted at least one infective, they are infected
      # l is the vector of people who were susceptible, who contacted at least
      # on infectious person
      ########################################################################
      l = which(vstate==0&vnum_infected_people_contacted>0)
      if (length(l)>0){
        vnewstate[l] = 1 

        ######################################################################
        # now, from the currently infectious people (those with vstate==1)
        # pick the one who infected each susceptible-but-now-infected person
        # zwho is the index of the person who infected the susceptible
        # ztime is the time at which that occurred
        # zlist is the list of secondary_infections of a particular infected person
        ######################################################################
        vprob = rep(0,length(vstate))
        vprob[vstate==1] = 1
        for (ll in l){
          m_parent = sample(seq(1,length(vstate)),1,prob=vprob,replace=T)
          zwho[ll] = m_parent
          ztime[ll] = t
          istage = 1
          while (length(m_parent)==1&m_parent>0){
            znum_descendant_infections[m_parent] = znum_descendant_infections[m_parent] + 1
            zlist[[m_parent]] = c(zlist[[m_parent]],ll)
            zgeneration[[m_parent]] = c(zgeneration[[m_parent]],istage)
            m_parent = zwho[m_parent]
            istage = istage + 1
          }
        }
      }

      vstate = vnewstate # update the state vector

      t = t + deltat     # update the time
   } 

   ###########################################################################
   # calculate the final size
   ###########################################################################
   final = 0
   if (length(vS)>0) final = max(vS)/N-min(vS)/N

   ###########################################################################
   # ztime_of_infection is the time at which a person was infected
   # zwho is the index of the person who infected them
   ###########################################################################

   l = which(z_linfected>0)
   wmax_generation = numeric(0)
   wnum_descendant_infections = numeric(0)
   for (i in l){
     amax = 0
     if (length(zlist[[i]])>0) amax = max(zgeneration[[i]])
     wmax_generation = c(wmax_generation,amax)
     wnum_descendant_infections = c(wnum_descendant_infections,length(zlist[[i]]))
   }

   return(list(time=vtime,I=vI,S=vS,final_size=final,who_infected_whom=zwho,time_of_infection=ztime,list_of_descendants=zlist,num_generations=zgeneration,linfected=z_linfected,num_descendant_infections=znum_descendant_infections,ANDI=mean(wnum_descendant_infections),avg_num_generations=mean(wmax_generation)))

} # end SIR_agent function definition


