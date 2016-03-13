number.of.each.age <- function(demographic.ages,disease.state,num.comps){
  
  number.of.age = matrix(0,length(demographic.ages[,1]),1)
  for (i in 1:length(demographic.ages[,1]))
  {
    age = demographic.ages[i,1]
    number.of.age[i]          =    sum(disease.state[((age*num.comps)+1):((age+1)*num.comps)])
  }
  return(number.of.age)
}

full.mixing.matrix<- function(contacts,demographic.ages){
  polymod.ages                =    contacts[,1]
  contacts                    <-   contacts[,2:16]
  max.age                     <-   max(demographic.ages[,1])
  contact.matrix              =    matrix(0,max.age+1,max.age+1)
  colnames(contact.matrix)    =    c(0:max.age)
  rownames(contact.matrix)    =    c(0:max.age)
  num.polymod.ages            <-   length(polymod.ages)
  row.start.position          =    1
  
  for(i in 1 : ( max(polymod.ages))){
    for (j in 1: (max(polymod.ages))){
      i1 = floor((i-1) / 5) + 1
      j1 = floor((j-1) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/5
    }
  }
  
  for(i in 1 : ( max(polymod.ages))){
    for ( j in (max(polymod.ages) + 1) : (max.age+1)){
      i1 = floor((i-1) / 5) + 1
      j1 = floor(max(polymod.ages) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/5
    }
  }
  
  for ( i in (max(polymod.ages) + 1) : (max.age+1)){
    for ( j in 1 :  ( max(polymod.ages))){
      i1 = floor(max(polymod.ages) / 5) + 1
      j1 = floor((j-1) / 5) + 1
      contact.matrix[ i, j]  =  contacts[ i1, j1]/ 23
    }
  }
  
  last.entry.col  =  length(contacts[1, ])
  last.entry.row  =  length(contacts[ , 1])
  for ( i in (max(polymod.ages) + 1) : (max.age+1)){
    for ( j in (max(polymod.ages) + 1) : (max.age+1)){
      contact.matrix[ i, j]  =  contacts[ last.entry.row, last.entry.col]/ 23
    }
  }
  
  
  return(contact.matrix)
}

average.contacts.by.age<- function(mixing.matrix,max.age){
  
  contacts.by.age       =     colSums(mixing.matrix[,seq(1,(max.age+1))])
  av.num.contacts       =     mean(contacts.by.age)
  return(av.num.contacts)
}
deterministic.transmission.matrix<-function(age , disease.state , vacc.immune , foi  , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*3)+1):((age+1)*3)]
  number.of.age             =    sum(age.disease.state)
  u                         =    time.step/365
  change.matrix             =    matrix(0,11,1)
  new.infecteds             =    0
  
  A                         =    matrix(0,11,1)
  A[1]                      =    round((1-u)*(1-foi) * age.disease.state[1])
  A[2]                      =    round((1-u)*(foi) * age.disease.state[1])   +   round((1-u)*(1-sigma) * age.disease.state[2])
  A[3]                      =    round((1-u)*(sigma) * age.disease.state[2])   +   round((1-u)*(1-rho) * age.disease.state[3])
  A[4]                      =    round((1-u)*(rho) * age.disease.state[3])   +   round((1-u) * age.disease.state[4])
  A[5]                      =    round((1-u)   *  age.disease.state[5])
  A[6]                      =    round((u)*(1-foi) * age.disease.state[1])
  A[7]                      =    round((u)*(foi) * age.disease.state[1])   +   round((u)*(1-sigma) * age.disease.state[2])
  A[8]                      =    round((u)*(sigma) * age.disease.state[2])   +   round((u)*(1-rho) * age.disease.state[3])
  A[9]                      =    round((u)*(rho) * age.disease.state[3])   +   round((u) * age.disease.state[4])
  A[10]                     =    round((u)   *  age.disease.state[5])
  
  A[11]                     =    round((1-u)*(foi) * age.disease.state[1])   +   round((u)*(foi) * age.disease.state[1])
  
  return(A)
}


stochastic.transmission.matrix<-function(age , disease.state , vacc.immune , foi , maternal.immunity.loss , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*num.comps)+1):((age+1)*num.comps)]
 # number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,2*num.comps,1)
  new.infecteds             =       0
  
  # Assume that no one 1 or older gets vaccinnated
  if (age == 1){
          v1         =   vacc.immune
  } else {
          v1         =    0
} 
          u          =    prob.age.change
          d          =    maternal.immunity.loss
  mat.immune         =    susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,12,1)

  if(age.disease.state[1] > 0)  {
    mat.immune       =    rmultinom(1 , age.disease.state[1] , c((1-u)*(1-d) , (1-u)*d*(1-v1) , 0 , 0 , 0 , (1-u)*d*v1 , 0 , u*(1-v1) , 0 , 0 , 0 ,u*v1))
  }else{
    mat.immune       =    matrix(0,2*num.comps,1)
  }
  
  if (age.disease.state[2] > 0){
      susceptibles     =    rmultinom(1, age.disease.state[2] , c(0, (1-u)*(1-foi)*(1-v1) , (1-u)*foi*(1-v1) , 0 , 0 , v1*(1-u) , 0 , u*(1-foi)*(1-v1) , u*foi*(1-v1) , 0 , 0 , u*v1))   
    
}
 
  if (age.disease.state[3] > 0){
    exposed     =    rmultinom(1, age.disease.state[3] , c(0, 0, (1-u)*(1-sigma) , (1-u)*sigma , 0 , 0 , 0 , 0 , u*(1-sigma) , u*sigma , 0 , 0))
    new.infecteds    =    exposed [4]    +    exposed[11]    ###### state S to state E 
    }
  
  if (age.disease.state[4] > 0){
    #number.recover   =    rbinom(age.disease.state[4],1,min(1,365*time.step/15))
    #prop.recover     =    number.recover/age.disease.state[4]
    infecteds        =    rmultinom(1,age.disease.state[4],c(0 , 0 , 0 , (1-u)*(1-rho) , (1-u)*rho , 0 , 0 , 0 , 0 , u*(1-rho) , u*rho,0))
  }
  
  if (age.disease.state[5] > 0){
    recovered        =    rmultinom(1,age.disease.state[5],c(0 , 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , 0 , u , 0))
  }
  
  if (age.disease.state[6] > 0){
    vac              =    rmultinom(1,age.disease.state[6],c(0 , 0 , 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , 0 , u))
  }
  
  change.matrix      =    c(mat.immune[1:6] + susceptibles[1:6] + exposed[1:6] + infecteds[1:6] + recovered[1:6] + vac[1:6],
                            mat.immune[7:12] + susceptibles[7:12] + exposed[7:12] + infecteds[7:12] + recovered[7:12] + vac[7:12])
  change.matrix[13] = new.infecteds
  return(change.matrix)
}

contacts.per.age.group <- function(mixing.matrix,demographic.ages){
  mixing.times.population <- mixing.matrix
  contacts.per.age <- mixing.matrix
  for (i in 1:length(demographic.ages[,1])){
    k<-matrix(demographic.ages[i,2],length(demographic.ages[,1]),1)
    mixing.times.population[,i] <- k*mixing.matrix[,i]
  }
  for (i in 1:length(demographic.ages[,1])){
    contacts.per.age[i,] <- mixing.times.population[i,]/demographic.ages[i,2]
  }
  
  return(rowSums(contacts.per.age))
}

grouped.ages <-function(contacts,demographic.ages)
{
  group.ages <- matrix(0,length(contacts[,1]),2)
  group.ages[,1] <- contacts[,1]
  count = 1
  for (i in 1:length(demographic.ages[,1]))
  {
    if(demographic.ages[i,1] > tail(group.ages[,1],1)){
      group.ages[length(group.ages[,1]),2] <- group.ages[length(group.ages[,1]),2] + demographic.ages[i,2]
    }
    else if(group.ages[count+1,1] > demographic.ages[i,1]){
      group.ages[count,2] <- group.ages[count,2] + demographic.ages[i,2] 
    }
    else {
      count <- count + 1
      group.ages[count,2] <- group.ages[count,2] + demographic.ages[i,2] 
    }
  }
  return(group.ages)
}

refresh.plots  <- function(j, infecteds.by.time, difference.from.estimate, foi.by.time, demographic.ages, prop.infected.by.age, prop.susceptible.age, t, pop.by.year, prop.sus.time, total.prop.susceptible ){
  par(mfrow=c(3,3))
  plot(1:j,infecteds.by.time[1:j],type="l",ylab = "new inf",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,difference.from.estimate[1:j],type="l",ylab = "diff.est",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,cumsum(infecteds.by.time[1:j]),type="l",ylab = "cum inf",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,total.prop.susceptible[1:j],type="l",ylab = "total prop sus",col=rgb(runif(1),runif(1),runif(1)))
  
  plot(demographic.ages[,1], prop.infected.by.age,type="l",ylab = "inf prop",xlab = "age",col=rgb(runif(1),runif(1),runif(1)))
  plot(1:j,prop.sus.time[1:j], type="b",ylab = "prop contacts sus",xlab = "step",col=rgb(runif(1),runif(1),runif(1)))
  plot(demographic.ages[,1],prop.susceptible.age,type="b",ylab = "prop susceptible",xlab = "age",col=rgb(runif(1),runif(1),runif(1)))
  years   =   ceiling(t/365)
  if (years > 0)
  {
    infections.per.year        =    matrix(0,years,1)
    for (pp in 1:years)
    {
      infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
    }
    plot(1:years, infections.per.year, type="b", ylab = "num infs", xlab = "year", col= rgb(runif(1), runif(1), runif(1)))
    plot(1:years, 100*infections.per.year/pop.by.year[1:years], type="b", ylab = "% of pop infected", xlab = "year", col= rgb(runif(1), runif(1), runif(1)))
  }
  
}

foi.by.next.gen <- function (mixing.matrix, disease.state, inf.comp, time.step , infectious.period, beta, demographic.ages, num.comps){
  number.age.brackets        =      length(disease.state)/num.comps
  infectious.indices         =      seq(inf.comp,length(disease.state),num.comps)
  number.infectious.by.age   =      disease.state[infectious.indices]
  pop.by.age  =  number.of.each.age (demographic.ages, disease.state, num.comps)
  foi.by.age  =  (colSums((number.infectious.by.age * t( mixing.matrix) * beta) ) / pop.by.age) * min( 1, time.step / infectious.period)
  return(foi.by.age)  
}

force.of.infection.by.age <- function(age,mixing.matrix,disease.state,beta,gamma,time.step,inf.comp,num.comps){
  number.age.brackets        =      length(disease.state)/num.comps
  infectious.indices         =      seq(inf.comp,length(disease.state),num.comps)
  number.infectious.by.age   =      disease.state[infectious.indices]
  mixing.by.age              =      mixing.matrix[,age+1]*time.step
  population.by.age          =      number.of.each.age (demographic.ages, disease.state, num.comps)
  foi.by.age                 =      1 - exp(-sum  ( beta* ( (number.infectious.by.age)^gamma )*mixing.by.age/population.by.age  )  )
  
  return(foi.by.age)
}

calibrate.beta <- function (mixing.matrix, disease.state, infectious.indices, max.age, time.step, infectious.period, R_0, population.by.age){
  num.infectious  =  matrix(0, (max.age + 1), 1)
  average.infectious  =  matrix(0, (max.age + 1), 1)
  for ( i in 1 : (max.age + 1) ){
    num.infectious[i]  =  1
    average.infectious[i]  =  sum(mixing.matrix[, i]) * min(1, time.step / infectious.period)
    num.infectious[i]  =  0
  }
  g = 0
  for (i in 1 : ( max.age + 1)){
    g = g + (population.by.age[i] / sum(population.by.age[1 : (max.age + 1)])) * average.infectious[i]
  }
  beta  = R_0 * min(1, time.step / infectious.period) / g
  return(beta)
}


initial.disease.state <- function(demographic.ages  ,  vacc.immune  ,  maternal.immunity.loss  ,  average.age.at.vaccination  ,  initial.prop.susceptible  ,  num.comps){
  disease.state                 =      matrix(0,num.comps*length(demographic.ages[,1]),1)
  disease.state[1]              =      ceiling(demographic.ages[1,2]*(1-maternal.immunity.loss))
  disease.state[2]              =      ceiling(demographic.ages[1,2]*(maternal.immunity.loss)*((average.age.at.vaccination)*(1-vacc.immune) + (1-average.age.at.vaccination)))
  disease.state[num.comps]      =      ceiling(demographic.ages[1,2]*(maternal.immunity.loss)*average.age.at.vaccination*(vacc.immune))
  
  for (i in 2:length(demographic.ages[,1])){
   disease.state[(((i-1)*num.comps)+1):(i*num.comps)]  =   ceiling(demographic.ages[i,2]*c(0,initial.prop.susceptible[i], initial.prop.exposed[i],
                     initial.prop.infected[i],initial.prop.recovery[i],
                     initial.prop.vac[i]))
}
  
  return(disease.state)
}


