setwd("/home/Mumps")

######################################################################################
######################################################################################
##########################   Full model: lag #########################################
######################################################################################
######################################################################################

library(pomp);require(magrittr)
load("dataset.RData")

#######################################################################################
################################# some structure for pomp #############################
#######################################################################################

rproc <- Csnippet("
                  double beta, br, seas, foi, gh, beta0, beta1, dw, births;
                  double rate[7], trans[7];
                  
                  // cohort effect
                  if (fabs(t-floor(t)-240.0/365.0) < 0.5*dt)
                  br = cohort*birthrate/dt + (1-cohort)*birthrate;
                  else
                  br = (1.0-cohort)*birthrate;
                  
                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  // wuyi: 5.1-5.3; shiyi: 10.1-10.7;summer:7.1-8.31;yuandan:1.1
                  if ( (t>=90&&t<=121) || (t>=274&&t<=304) || (t>=4&&t<40) || (t>=124 && t<212) || (t>= 314 && t<365))
                  seas = 1.0+amplitude*0.34/0.66;
                  else
                  seas = 1.0-amplitude;
                  
                  // transmission rate
                  beta0 = R0*gamma*seas;
                  beta1 = exp(bmt*MT);
                  //beta1 = pow(2,bmt*MT); 
                  beta = beta0*beta1;        
                  
                  // expected force of infection
                  gh  = pow(I+iota,alpha);
                  foi = beta*gh/pop;    // add the effect of temp on suspect and
                  
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  rate[0] = foi*dw/dt;      // stochastic force of infection
                  rate[1] = mu;             // natural S death
                  rate[2] = sigma;        // rate of ending of latent stage
                  rate[3] = mu;             // natural E death
                  rate[4] = gamma;          // recovery
                  rate[5] = mu;             // natural I death
                  rate[6] = mu;             // natural R death
                  
                  // Poisson births
                  births = rpois(br*dt*0.22);  // 0.7= 0.984 * 0.8
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  reulermultinom(1,R,&rate[6],dt,&trans[6]);
                  
                  S += births   - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  if (!R_FINITE(S)) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg \\\n\",R0,beta0,beta1,beta,foi,I,alpha);
                  ")



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
                  double rhonew = exp(MT*gmt)*rho/(exp(MT*gmt)*rho+1);
                  double m = rhonew * C;
                  double v = m*(1.0-rhonew+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
                  double rhonew = exp(MT*gmt)*rho/(exp(MT*gmt)*rho+1);
                  double m = rhonew*C;
                  double v = m*(1.0-rhonew+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

toEst <- Csnippet("
                  Tmu = logit(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tcohort = logit(cohort);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  TR0 = log(R0);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = expit(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    Tcohort = expit(cohort);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    TR0 = exp(R0);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")

#######################################################################################
################################# Gather pomp structure ###############################
#######################################################################################

mumps.hz %>%
  pomp(t0=mumps.hz$week[1],
       time="week",
       rprocess=discrete.time.sim(step.fun=rproc,delta.t=3/365.25),
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covar,
       toEstimationScale=toEst,
       fromEstimationScale=fromEst ,
       tcovar="week",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "bmt","gmt",
                    "S_0","E_0","I_0","R_0")
  ) -> origin.model

#######################################################################################
#################################### Simulation #######################################
#######################################################################################

require(foreach);require(doMC);
registerDoMC(16)
set.seed(998468235L,kind="L'Ecuyer")
mcopts <- list(preschedule=FALSE,set.seed=TRUE)


theta.lo <- c(mu=5.93/1000,sigma=5,gamma=5,alpha=3.115105e-01,iota=5.952648e+00,rho=0.2,R0=20,
              sigmaSE=0,psi=0.1,cohort=0,amplitude=0, S_0=0.02,I_0=2e-06, E_0=2e-06,bmt=-2,gmt=-2)
theta.hi <- c(mu=5.93/1000,sigma=40,gamma=40,alpha=1,iota=40,rho=0.7,R0=200,
              sigmaSE=1,psi=0.4,cohort=1,amplitude=1, S_0=0.03,I_0=0.0001,E_0=0.0001,bmt=2,gmt=2)
gusses.parameter <-cbind(theta.lo,theta.hi)

rw.sd_rp <- 0.02 ; rw.sd_cli <- 0.03; rw.sd_ivp <-0.02


stew(file=sprintf("globalltemp1temp3.rda"),{
  global.result<- foreach(i=1:2000, .combine=rbind,
                          .packages=c("pomp","magrittr"),.errorhandling="remove",
                          .inorder=FALSE, .options.multicore=mcopts) %dopar%  {
                            start.parameter <- apply(gusses.parameter,1,function(x)runif(1,x[1],x[2]))
                            start.parameter[["R_0"]]<-1-start.parameter[["S_0"]]-start.parameter[["E_0"]]-start.parameter[["I_0"]]
                            mif2(origin.model, start=start.parameter, Np=2000,Nmif=20,
                                 cooling.type="geometric",cooling.fraction.50=0.1,transform=TRUE,
                                 rw.sd=rw.sd(
                                   gamma=rw.sd_rp,alpha=rw.sd_rp,iota=rw.sd_rp, sigma=rw.sd_rp,
                                   sigmaSE=rw.sd_rp, psi=rw.sd_rp, cohort=rw.sd_rp,R0=rw.sd_rp,
                                   amplitude=rw.sd_rp,gamma=rw.sd_rp,
                                   bmt=rw.sd_cli,gmt=rw.sd_cli,
                                   I_0=ivp(rw.sd_ivp),
                                   S_0=ivp(rw.sd_ivp),
                                   E_0=ivp(rw.sd_ivp),
                                   R_0=ivp(rw.sd_ivp)
                                 )) -> mf   ##  %>%  mif2()
                            pf <- replicate(10, pfilter(mf, Np =2000))
                            ll <- sapply(pf,logLik)
                            ll <- logmeanexp(ll, se = TRUE)
                            nfail <- sapply(pf,getElement,"nfail")
                            data.frame(as.list(coef(mf)),
                                       loglik = ll[1],
                                       loglik.se = ll[2],
                                       nfail.min = min(nfail),
                                       nfail.max = max(nfail))
                          }},seed=290860873,kind="L'Ecuyer")
write.csv(global.result,file="globaltemp1temp3.csv",row.names=F)   ### true

