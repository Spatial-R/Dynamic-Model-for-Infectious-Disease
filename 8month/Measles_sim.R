#library("prob")
#library("ggplot2")
#library("RColorBrewer")
source("Measles_func.R")

demographic.ages              =       read.csv("age1.csv")
contacts                      <-      read.csv("Contacts.csv")
intial                        <-      read.csv("initial.csv")

num.comps                     =       6          # number of compartments in the model
inf.comp                      =       4          # what compartment keeps track of the number of infectious for each age group
t                             =       0
R_0                           =       4.2
time.step                     =       5        # number of days to step forward
average.age.at.vaccination    =       0.67
vacc.prop                     =       0.90
vacc.success                  =       0.80
vacc.immune                   =       vacc.prop*vacc.success*(1/12)
#v                             =       vacc.success*min(1,time.step/60)*min(11/12,(1+time.step/60)/6)
maternal.immunity.loss        =       min(1,time.step/250)
birth.rate                    =       10.51/(1000*365)
death.rate                    =       5.51/(1000*365)
prob.survival                 =       1-(death.rate)*(time.step)
initial.prop.susceptible      =       intial[,3]*2    ##### the percentage of susceptible
initial.prop.exposed          =       intial[,4]  ##### the percentage of exposed
initial.prop.infected         =       intial[,5]  ##### the percentage of infected
initial.prop.vac              =       intial[,7]+intial[,8]     ##### the percentage of vac1
initial.prop.recovery         =   1-intial[,2]-intial[,3]-intial[,4]-intial[,5]-intial[,7] 

years.per.extra.vaccination   =       3
gamma                         =       0.97
exposed.days                  =       18       # number of days spent in the exposed class on average  14~18
infectious.days               =       8       # number of days spent in the infected class on average 6-8
sigma                         =       min(1,time.step/exposed.days)     ## the rate at which individuals move from the exposed to the infectious classes. Its reciprocal (1/σ) is the average latent (exposed) period.
rho                           =       min(1,time.step/infectious.days)  ## its reciprocal (1/γ) which determines the average infectious period.
max.age                       =       40      # assume that when R_0 was calculated previously, this was the maximum age of the people spreading measles

mixing.matrix                 <-      full.mixing.matrix(contacts,demographic.ages)      # average number of people met of each age grooup, stratified by age
mean.number.contacts          <-      average.contacts.by.age(mixing.matrix,max.age)
beta_0                        =       R_0/(mean.number.contacts*infectious.days)         # average number of people infected during infectiousness = R_0, 
                                                                                         # therefore mean transmission rate is given by this expression
beta_1                        =       0.8

disease.state                 <-      initial.disease.state(demographic.ages,  vacc.immune  ,  maternal.immunity.loss  ,  average.age.at.vaccination  ,  initial.prop.susceptible  ,  num.comps)

av.migrants.per.age.per.day   =       1/(365 * length(demographic.ages[,1])) 
prob.survival                 =       1 - death.rate*time.step
updated.state                 =       matrix(0,num.comps*length(demographic.ages[,1]),1)
num.steps                     =       4000
infecteds.by.time             =       matrix(0,num.steps,1)
susceptibles.by.time          =       matrix(0,num.steps,1)
average.infection.age         =       matrix(0,num.steps,1)

# make transition matrix for all ages
for(j in 1:num.steps){
  N                           <-      sum(disease.state)
  new.infected                =       0
  beta                        =       beta_0 * ( 1 + beta_1 * cos( 2 * pi * (j*time.step)) )
  average.births              =       birth.rate*N*time.step
  total.births                =       rpois(1,average.births)
  disease.state[1]            =       disease.state[1] + total.births
  migrant.infecteds           =       rpois(length(demographic.ages[,1]),av.migrants.per.age.per.day*time.step)
  #print(paste("infected.migrants =",sum(migrant.infecteds)))
  disease.state[seq(inf.comp,length(disease.state),num.comps)]    =  disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )] + migrant.infecteds
  #print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
  
 updated.state               =       matrix( 0, num.comps*length(demographic.ages[,1]), 1)

  for ( i in 1:length(demographic.ages[,1]) ){
     age                     =      demographic.ages[i,1]
     foi.age                 =      force.of.infection.by.age(age,mixing.matrix,disease.state,beta,gamma,time.step,inf.comp,num.comps)
     change.matrix.by.age    =      stochastic.transmission.matrix(age , disease.state , vacc.immune , foi.age , maternal.immunity.loss , demographic.ages , time.step  ,sigma  ,  rho)
     new.infected            =      new.infected   +    change.matrix.by.age[(2*num.comps+1)]
   if(age == max(demographic.ages[,1])){
        updated.state[seq((((i-1)*num.comps)+1),num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      } else{
        updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[1:(2*num.comps)]
      }
    }

  total.infecteds                =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)])
  number.by.age                  =       number.of.each.age(demographic.ages,updated.state,num.comps)
  prop.infected.by.age           =       updated.state[seq(inf.comp,length(updated.state),num.comps)]/number.by.age
  average.infection.age[j]       =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)]*demographic.ages[,1])
  infecteds.by.time[j]           =       new.infected
  susceptibles.by.time[j]        =       sum(updated.state[seq(2,length(updated.state),num.comps)])
  disease.state                  =       ceiling(updated.state*prob.survival) 
        t                        =       t   +    time.step
     
        if (j %% 50 == 0){
          print(paste("j =",j))
          # Sys.sleep(0.01)
          par(mfrow=c(2,3))
          plot(1:j,infecteds.by.time[1:j],type="l",ylab = "new inf",xlab="step")
          plot(1:j,cumsum(infecteds.by.time[1:j]),type="l",ylab = "cum inf",xlab ="step")
          plot(demographic.ages[,1],prop.infected.by.age,type="l",ylab = "inf prop",xlab = "age")
        #  plot(1:j,all.infectious.by.time[1:j],type="l",ylab = "all inf")
          years   =   ceiling(t/365)
          if (years > 0) {
            infections.per.year        =    matrix(0,years,1)
            for(pp in 1:years){
              infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
            }
            plot(1:years,infections.per.year,type="b",ylab = "num infs",xlab = "year")
          }
        }
}

matrix(disease.state,ncol=6,byrow = T)
