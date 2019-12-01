#analysis.tem.3 <- analysis.tem.3[,c(2,5:10,17)]
#analysis.tem.3$region <- factor(analysis.tem.3$region,levels =  unique(analysis.tem.3$region),
#                                labels = c(letters[1:26],LETTERS[1:3]))

library(brms)
setwd("E:/Github/Dynamic-Model-for-Infectious-Disease/Example")
load("test.RData")

prior <- get_prior(bf(log(adjust.rc) ~ s(AH) + s(rain) + s(sunshine) + (1 | season) + (1 | region)),
          data = analysis.tem.3)
prior$prior[1] <- prior$prior[5] <- prior$prior[6] <- prior$prior[11] <- prior$prior[15] <- "student_t(10,0,1)"
 
fit6 <- brm(bf(log(adjust.rc) ~ s(AH) + s(rain) + s(sunshine) + (1 | season) + (1 | region)),
            data = analysis.tem.3,prior = prior,chains = 6,cores = 6,
            control = list(adapt_delta = 0.99),iter = 2000)
marginal_smooths(fit6)


