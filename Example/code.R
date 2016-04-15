setwd("/home/varicella")
library(pomp)
library(magrittr)
load("dataset.RData")

rproc <- Csnippet("
                  double Beta, lambda, dW, airplt, logbeta1;
                  double rate[7], trans[7];
                
                  logbeta1 = log(R0*(1.0-exp(-(gamma+mu)*dt))/dt);    // caculate the log beta based on R0

                  //  seasonality for transmission
                  int nseas=3;
                  int k;
                  const double *seas = &seas1;
                  const double *logbeta = &logbeta1;
                  for (k = 0, Beta = 0; k < nseas; k++) Beta += logbeta[k] *seas[k];

                  //  air polutant on  transmission
                  airplt =  pm10 * pam1pm10+ pam1temp * temp ;
                  Beta =exp(Beta * airplt+ effcontact);    // add the contact infection

                  dW = rgammawn(sigmaSE,dt);   // white noise (extrademographic stochasticity)

                  // expected force of infection
                  lambda = (Beta*pow(I+iota,alpha)/pop)*(dW/dt);   // mixing pattern and immigration 
                  

                  double   birthnum    = rpois(birthrate*dt);		     // births are Poisson
                  double   vaccinum    = rpois(vacci*dt);      // vaccination are Poisson

                  rate[0] = (pm10 * pam2pm10+ pam2temp * temp)*lambda * dt;  // stochastic force of infection with the effects on susceptible
                  rate[1] = mu * dt;             // natural S death
                  rate[2] = sigma * dt;        // rate of ending of latent stage
                  rate[3] = mu * dt;             // natural E death
                  rate[4] = (pm10 * pam3pm10+ pam3temp * temp)*gamma * dt;        // recovery
                  rate[5] = mu * dt;             // natural I death
                  rate[6] = mu * dt;             // natural R death
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  reulermultinom(1,R,&rate[6],dt,&trans[6]);
                  
                  S += birthnum - trans[0] - trans[1] - vaccinum ;
                  E += trans[0]  - trans[2] - trans[3];
                  I += trans[2]  - trans[4] - trans[5];
                  R += trans[4]  - trans[6] + vaccinum; //pop-S-E-I;
                  W += (dW- dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  ")


Csnippet("
                 double Beta, lambda, airplt, logbeta1;
                 double rate[7], trans[7];
         
                 logbeta1 = log(R0*(1.0-exp(-(gamma+mu))));    // caculate the log beta based on R0
         
                  //  seasonality for transmission
                 int nseas=3, k;
                 const double *seas = &seas1;
                 const double *logbeta = &logbeta1;
                 for (k = 0, Beta = 0; k < nseas; k++) Beta += logbeta[k] *seas[k];
         
                //  air polutant on  transmission
                airplt = pm10 * pam1pm10+ pam1temp * temp;
                Beta =exp(Beta * airplt + effcontact);
         
                // expected force of infection
                lambda = Beta*pow(I+iota,alpha)/pop;   // mixing pattern and immigration 
         
                double  birthnum    = birthrate;		     // births are Poisson
                double  vaccinum    = vacci;      // vaccination are Poisson
         
                rate[0] =  (pm10 * pam2pm10+ pam2temp * temp)*lambda  ;  // stochastic force of infection
                rate[1] = mu ;             // natural S death
                rate[2] = sigma  ;        // rate of ending of latent stage
                rate[3] = mu  ;             // natural E death
                rate[4] = (pm10 * pam3pm10+ pam3temp * temp)*gamma  ;        // recovery
                rate[5] = mu ;             // natural I death
                rate[6] = mu ;             // natural R death
         
               // transitions between classes

               trans[0] = S * rate[0];
               trans[1] = S * rate[1];
               trans[2] = E * rate[2];
               trans[3] = E * rate[3];
               trans[4] = I * rate[4];
               trans[5] = I * rate[5];
               trans[6] = R * rate[6];
         
               DS = birthnum  - trans[0] - trans[1] - vaccinum ;
               DE = trans[0]  - trans[2] - trans[3];
               DI = trans[2]  - trans[4] - trans[5];
               DR = trans[4]  - trans[6] + vaccinum;
               DC = trans[4];           // true incidence
         ") -> skel

initlz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")

dmeas <- Csnippet("
                  double m =  rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tpsi = log(psi);
                  Teffcontact= log(effcontact);
                  //Tpamso2 = log(pamso2);
                  //Tpamno2 = log(pamno2);
                  //Tpampm10 = log(pampm10);
                  //Tpampm25 = log(pampm25);
                  //Tpamco = log(pamco);
                  //Tpamo3 = log(pamo3);
                  TsigmaSE = log(sigmaSE);
                  TR0 = log(R0);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tpsi = exp(psi);
                    Tiota = exp(iota);
                    Teffcontact=exp(effcontact);
                    Trho = expit(rho);
                    //Tpamso2 = exp(pamso2);
                    //Tpamno2 = exp(pamno2);
                    //Tpampm10 = exp(pampm10);
                    //Tpampm25 = exp(pampm25);
                    //Tpamco = exp(pamco);
                    //Tpamo3 = exp(pamo3);
                    TsigmaSE = exp(sigmaSE);
                    TR0 = exp(R0);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")


c(R0=5,gamma=1/20,sigma=1/15,mu=1.28e-04,sigmaSE=0.512,
  rho=0.8, alpha= 0.976, iota=2.9, psi=0.116,effcontact=0.5,
  pam1pm10=0.0001,pam1temp=0.002, #pam1pm25=0.002
  pam2pm10=0.0001,pam2temp=0.002, #pam2pm25=0.002,
  pam3pm10=0.0001,pam3temp=0.002, #pam3pm25=0.002,
  S_0=0.0297,I_0=5.14e-05,E_0=5.17e-05
) -> theta
theta["R_0"]<-1-theta["I_0"]-theta["E_0"]-theta["E_0"]

#######################################################################################
################################# Gather pomp structure ###############################
#######################################################################################

varicella.hz %>%
  pomp(t0=varicella.hz$week[1],
       time="week",
       params=theta,
       rprocess=discrete.time.sim(step.fun=rproc,delta.t=1),
       skeleton=skel,skeleton.type="map", skelmap.delta.t=1,
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covar,
       toEstimationScale=toEst,
       fromEstimationScale=fromEst ,
       tcovar="week",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota","psi","effcontact",
                    "rho","sigmaSE","pam1pm10","pam1temp",  ###,"pam1pm25",
                    "pam2pm10","pam2temp",                  ### "pam2pm25",
                    "pam3pm10","pam3temp",                  ### "pam3pm25",
                    "S_0","E_0","I_0","R_0")
  ) -> origin.model



require(foreach);require(doMC); 
registerDoMC(4)
set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)

run_level <- 2
mumps_Np <-          c(100,5e3,1e4)
mumps_Nmif <-        c(10, 100,400)
mumps_Nreps_eval <-  c(2,  5,  20)
mumps_Nreps_local <- c(10, 20,  40)
mumps_Nreps_global <-c(10, 20, 100)
mumps_Nsim <-        c(50,100, 500) 


gusses.parameter <- rbind( 
  R0=c(0,30),
  mu=c(0.002,0.020),
  sigma=c(0,25),
  gamma=c(0,35),
  alpha=c(0.976,0.976),
  iota=c(0,10),
  rho=c(0,1),
  sigmaSE=c(0,1),
  psi=c(0,1),
  pam1pm10=c(-0.0001,0.0001),
  pam1temp=c(-0.001,0.001),
  pam2pm10=c(-0.0001,0.0001),
  pam2temp=c(-0.001,0.001),
  pam3pm10=c(-0.0001,0.0001),
  pam3temp=c(-0.001,0.001),
  effcontact=c(0,1),
  S_0=c(0,0.03),
  I_0=c(0,0.0001),
  E_0=c(0,0.0001)
)

stew(file=sprintf("global0617.rda",run_level),{
  global.result<- foreach(i=1:100, .combine=rbind,
                          .packages=c("pomp","magrittr"),.errorhandling="pass",
                          .inorder=FALSE, .options.multicore=mcopts) %dopar%  {
                            start.parameter<-apply(gusses.parameter,1,function(x)runif(1,x[1],x[2]))
                            start.parameter[["R_0"]]<-1-start.parameter[["S_0"]]-start.parameter[["E_0"]]-start.parameter[["I_0"]]
                            mif2(origin.model, start=start.parameter, Np=1000,
                                 Nmif=100,
                                 cooling.type="geometric",cooling.fraction.50=0.1, transform=TRUE,
                                 rw.sd=rw.sd( 
                                   R0=0.02,mu=0.02,sigma=0.02,gamma=0.02,
                                   alpha=0.02,iota=0.02, rho=0.02,
                                   sigmaSE=0.02, psi=0.02, pam1pm10=0.02,
                                   pam2pm10=0.02,pam3pm10=0.02, effcontact=0.02,
                                   pam1temp=0.02,
                                   pam2temp=0.02,
                                   pam2temp=0.02,
                                   I_0=ivp(0.02),
                                   S_0=ivp(0.02),
                                   E_0=ivp(0.02),
                                   R_0=ivp(0.02)
                                 )) -> mf   ##  %>%  mif2() 
                            pf <- replicate(5, pfilter(mf, Np =1000))
                            ll <- sapply(pf,logLik)
                            ll <- logmeanexp(ll, se = TRUE)
                            nfail <- sapply(pf,getElement,"nfail")
                            data.frame(as.list(coef(mf)),
                                       loglik = ll[1],
                                       loglik.se = ll[2],
                                       nfail.min = min(nfail),
                                       nfail.max = max(nfail))
                          }},seed=290860873,kind="L'Ecuyer")

write.csv(global.result,file="global.result.csv",row.names=F)
