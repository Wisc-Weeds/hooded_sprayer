set.seed(271828)
data(sleepstudy)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
fm1

sleepstudy %>% 
  distinct(Subject)

library(merTools)
PI.time <- system.time(
  PI <- predictInterval(merMod = fm1, newdata = sleepstudy,
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
)


library(ggplot2)
ggplot(aes(x=1:181, y=fit, ymin=lwr, ymax=upr), data=PI[1:181,]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") + theme_bw()


fm1 <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
display(fm1)


sleepstudy[1,]
#>   Reaction Days Subject
#> 1   249.56    0     308
predictInterval(fm1, sleepstudy[1,], include.resid.var=0) #predict the average body fat for a group of 196cm female baseball players


predictInterval(fm1, sleepstudy[1,], include.resid.var=0, ignore.fixed.terms = 1)

predictInterval(fm1, sleepstudy[1,], include.resid.var=0, ignore.fixed.terms = "(Intercept)")

predictInterval(fm1, sleepstudy[1,], include.resid.var=0, ignore.fixed.terms = 1:2)


predictInterval(fm1, sleepstudy[1,], include.resid.var=0,
                fix.intercept.variance = TRUE)




#Step 2: Comparison with arm::sim()



PI.arm.time <- system.time(
  PI.arm.sims <- arm::sim(fm1, 1000)
)

PI.arm <- data.frame(
  fit=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.500)),
  upr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.975)),
  lwr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.025))
)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="arm::sim()", x=(1:nrow(PI.arm))+0.1, PI.arm))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)




#Step 3a: lme4::bootMer() method 1


##Functions for bootMer() and objects
####Return predicted values from bootstrap
mySumm <- function(.) {
  predict(., newdata=sleepstudy, re.form=NULL)
}
####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

##lme4::bootMer() method 1
PI.boot1.time <- system.time(
  boot1 <- lme4::bootMer(fm1, mySumm, nsim=250, use.u=FALSE, type="parametric")
)

PI.boot1 <- sumBoot(boot1)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 1", x=(1:nrow(PI.boot1))+0.1, PI.boot1))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)


#####Step 3b: lme4::bootMer() method 2

##lme4::bootMer() method 2
PI.boot2.time <- system.time(
  boot2 <- lme4::bootMer(fm1, mySumm, nsim=250, use.u=TRUE, type="parametric")
)

PI.boot2 <- sumBoot(boot2)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 2", x=(1:nrow(PI.boot2))+0.1, PI.boot2))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)





# Step 3c: lme4::bootMer() method 3
##lme4::bootMer() method 3
PI.boot3.time <- system.time(
  boot3 <- lme4::bootMer(fm1, mySumm, nsim=250, use.u=TRUE, type="semiparametric")
)

PI.boot3 <- sumBoot(boot3)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 3", x=(1:nrow(PI.boot3))+0.1, PI.boot3))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)








#Step 3c: Comparison to rstanarm

library(rstanarm)
library(rstantools)
library(rstan)

PI.time.stan <- system.time({
  fm_stan <- stan_lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy,
                       verbose = FALSE, open_progress = FALSE, refresh = -1,
                       show_messages=FALSE, chains = 1)
  zed <- posterior_predict(fm_stan)
  PI.stan <- cbind(apply(zed, 2, median), predictive_interval(zed, prob=0.95))
})
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: 
#> Chain 1:  Elapsed Time: 6.994 seconds (Warm-up)
#> Chain 1:                2.497 seconds (Sampling)
#> Chain 1:                9.491 seconds (Total)
#> Chain 1:


print(fm_stan)
#> stan_lmer
#>  family:       gaussian [identity]
#>  formula:      Reaction ~ Days + (Days | Subject)
#>  observations: 180
#> ------
#>             Median MAD_SD
#> (Intercept) 251.5    6.4 
#> Days         10.5    1.7 
#> 
#> Auxiliary parameter(s):
#>       Median MAD_SD
#> sigma 25.9    1.6  
#> 
#> Error terms:
#>  Groups   Name        Std.Dev. Corr
#>  Subject  (Intercept) 23.8         
#>           Days         6.9     0.09
#>  Residual             26.0         
#> Num. levels: Subject 18 
#> 
#> ------
#> * For help interpreting the printed output see ?print.stanreg
#> * For info on the priors used see ?prior_summary.stanreg

PI.stan <- as.data.frame(PI.stan)
names(PI.stan) <- c("fit", "lwr", "upr")
PI.stan <- PI.stan[, c("fit", "upr", "lwr")]
comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="rstanArm", x=(1:nrow(PI.stan))+0.1, PI.stan))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)






fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
(p1 <- plotFEsim(FEsim(fm1)))


fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
(p1 <- plotREsim(REsim(fm1)))
#Plot just the random effects for the Days slope
(p2 <- plotREsim(REsim(fm1), facet= list(groupFctr= "Subject", term= "Days")))




fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:181,])
head(mfx)

ggplot(mfx) + aes(x = breaks, y = fit_Subject, group = case) +
  geom_line() +
  facet_wrap(~term)







InstEval
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data








m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
regFit <- predict(m1, newdata = sleepstudy[11, ]) # a single value is returned
intFit <- predictInterval(m1, newdata = sleepstudy[11, ]) # bounded values
# Can do glmer
d1 <- cbpp
d1$y <- d1$incidence / d1$size
gm2 <- glmer(y ~ period + (1 | herd), family = binomial, data = d1,
             nAGQ = 9, weights = d1$size)
regFit <- predict(gm2, newdata = d1[1:10, ])
# get probabilities
regFit <- predict(gm2, newdata = d1[1:10, ], type = "response")
intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "probability")
intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "linear.prediction")










#For a one-level random intercept model
m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
m1.er <- REimpact(m1, newdata = sleepstudy[1, ], breaks = 2)
#For a one-level random intercept model with multiple random terms
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#ranked by the random slope on Days
m2.er1 <- REimpact(m2, newdata = sleepstudy[1, ],
                   groupFctr = "Subject", term="Days")
#ranked by the random intercept
m2.er2 <- REimpact(m2, newdata = sleepstudy[1, ],
                   groupFctr = "Subject", term="int")
# You can also pass additional arguments to predictInterval through REimpact
g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
zed <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", n.sims = 50,
                include.resid.var = TRUE)
zed2 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "s", n.sims = 50,
                 include.resid.var = TRUE)
zed3 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", breaks = 5,
                 n.sims = 50, include.resid.var = TRUE)









library(merTools)
library(lme4)

dat <- iris
mod <- lmer(Sepal.Length ~ 1 + (1 + Sepal.Width + Petal.Length + 
                                  Petal.Width|Species), data=dat)

c1 <- predict(mod, dat)
c2 <- predictInterval(mod, dat)

plot_data <- cbind(c1, c2)
plot_data$order <- c(1:nrow(plot_data))

library(ggplot2)
ggplot(plot_data) + geom_line(aes(x=order, y=c1), color='red') + 
  geom_ribbon(aes(x=order, ymin=lwr, ymax=upr), color='blue', alpha=0.2) +  
  geom_line(aes(x=order, y=fit), color='blue')
