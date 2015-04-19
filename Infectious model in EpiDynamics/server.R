library(shiny);library(reshape)
library(ggplot2)
library(gridExtra)
library(scales)
library(EpiDynamics)

shinyServer(function(input, output, session) {  
  
  data<-reactive({

            if (input$modtype == "SIR"){  
    SIR(pars=c(beta=input$beta,gamma=input$gamma),
            init =c(S=input$S,I=input$I,R=input$R),time =0:input$time)$results
    } else if (input$modtype == "SIR + 2 Age Class") {
        sir2AgeClasses(
            pars=c(betaCC=input$betaCC,betaCA=input$betaCA,betaAC=input$betaAC,
              betaAA=input$betaAA,gamma=input$gamma,
                  lC=input$lC,muC=input$muC,muA=input$muA),
           init=c(SC=input$SC,IC=input$IC,SA=input$SA,IA=input$IA),time =c(1:input$time))$results
    } else if (input$modtype   ==  "SIR + 2 Type Imports") {
      SIR2TypesImports(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu,
               epsilon=input$epsilon,delta=input$delta),
        init=c(X=input$X,Y=input$Y,N=input$N),time=input$time)$results
    } else if (input$modtype   ==  "SIR + Constant Additive"){
      SIRAdditiveNoise(
        pars= c(beta=input$beta,gamma=input$gamma,mu=input$mu,
                noise=input$noise,N=input$N),
        init=c(X=input$X,Y=input$Y),time=c(0:input$time))$results
    } else if (input$modtype   ==  "SIR + Birth Death"){
      SIRBirthDeath(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu),
        init=c(S=input$S,I=input$I,R=(1-input$S-input$I)),time=c(0:input$time))$results
    } else if (input$modtype == "SIR + Carrier State"){
      SIRCarrierState(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu,
               epsilon=input$epsilon,rho=input$rho,Gamma=input$Gamma),
        init=c(S=input$S,I=input$I,C=input$C,R=(1-input$S-input$I-input$C)),
        time=c(0:input$time))$results
    } else if (input$modtype   ==  "SIR + DemogStoch") {
      SIRDemogStoch(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu),
        init=c(X=input$X,Y=input$Y,N=input$N),time=input$time)$results
    } else if (input$modtype == "SIR + Induced Mortality") {
      SIRInducedMortality2(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu,
               rho=input$rho,nu=input$nu),
        init= c(X=input$X,Y=input$Y,Z=input$Z),time=c(0:input$time))$results
    } else if (input$modtype == "SIR + Partial Immunity"){
      SIRPartialImmunity(
        pars=c(beta1=input$beta1,beta2=input$beta2,gamma1=input$gamma1,
               gamma2=input$gamma2,alpha1=input$alpha1,alpha2=input$alpha2,
               a1=input$a1,a2=input$a2,mu=input$mu,v=input$v),
        init=c(NSS=input$NSS,NIS=input$NIS,NRS=input$NRS,NRI=input$NRI,
               NSI=input$NSI,NSR=input$NSR,NIR=input$NIR,NRR=input$NRR),
        time=c(0:input$time))$results
    } else if  (input$modtype=="SIR + Sinusoidal Births"){
      SIRSinusoidalBirth(
        pars=list(beta=input$beta,gamma=input$gamma,alpha0=input$alpha0,
               mu=input$mu,alpha1=input$alpha1,w=input$w),
        init=c(S=input$S,I=input$I,R=1-(input$S+input$R)),time=c(0:input$time))$results
    } else if (input$modtype=="SIR + Sinusoidal Forcing"){
      SIRSinusoidalForcing(
        pars=list(beta1=input$beta1,beta0=input$beta0,gamma=input$gamma,
               omega=input$omega,mu=input$mu),
        init=c(S=input$S,I=input$I,R=1-(input$I+input$S)),time=0:200)$results
    } else if (input$modtype=="SIR + Tau Leap"){
      SIRTauLeap(
        pars=c(beta=input$beta,gamma=input$gamma,mu=input$mu,N=input$N,tau=input$tau),
        init=c(X=input$X,Y=input$Y,Z=input$Z),end.time=input$time)$results
    } else  {
      SIRVector(
        pars=c(muM=input$muM,muH=input$muH,vH=input$vH,vM=input$vM,betaHM=input$betaHM,
               betaMH=input$betaMH,gamma=input$gamma,r=input$r),
        init=c(XH=input$XH,XM=input$XM,YH=input$YH,YM=input$YM),time=c(0:input$time))$results
    }
  })    
  
  output$plot <- renderPlot({
     data.plot<-melt(data(),id.vars="time")
   ggplot(data.plot,aes(x=time,y=value,colour=variable)) + geom_line()+
       theme_bw(base_size=12,base_family="serif")+
      theme(legend.position = "top", legend.title = element_blank(),
            plot.margin = unit(c(0, 1, 1, 1), "line")) + 
      xlab("Time") + ylab("Number of people")   
  })
  output$table <- renderDataTable({
    k<-data();k[,2:dim(k)[2]]<-apply(k[,2:dim(k)[2]],2,round,4);k
  })
})


