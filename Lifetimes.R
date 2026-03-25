############################################################################################
#     
#     Set seed
#
RNGkind(sample.kind = "Rejection")
set.seed(5828) # change this to the last 4 digits of your student ID
#
#     Generate American Shorthair lifetimes with censoring
#
LenA <- sample(160:180, 1)
LifetimeA <- c(sort(round(pmin(rgamma(LenA-10,1.2,0.4), 15), 6)), rep(15,10))
CensorA <- c(rbinom(LenA,1,0.8))
CensorA[LifetimeA == 15] <- 0
#
#     Generate British Shorthair lifetimes with censoring
#
LenB <- sample(140:160, 1)
LifetimeB <- c(sort(round(pmin(rgamma(LenB-10,1,0.3), 15), 6)), rep(15,10))
CensorB <- c(rbinom(LenB,1,0.8))
CensorB[LifetimeB == 15] <- 0
#
#     Combine data & Generate ages at FIP contraction
#
Lifetime <- c(LifetimeA, LifetimeB)
Censor <- c(CensorA, CensorB)
Group <- factor(c(rep(1,LenA), rep(2,LenB)))
FIPAge <- ifelse(Lifetime>2, ifelse(runif(LenA+LenB)<0.4, pmax(Lifetime-rexp(LenA+LenB,1.5),0), 0), 0)
#
############################################################################################
# Please insert your R code after this line
############################################################################################
library(survival)

#Separate breeds
Adata = data.frame(time = Lifetime[Group==1],
                    status = Censor[Group==1])

Bdata = data.frame(time = Lifetime[Group==2],
                    status = Censor[Group==2])

Alldata = data.frame(time = Lifetime,
                      status = Censor)

fitA = survfit(Surv(time, status) ~ 1, data = Adata, conf.type="plain")
fitB = survfit(Surv(time, status) ~ 1, data = Bdata, conf.type="plain")
fitAll = survfit(Surv(time, status) ~ 1, data = Alldata, conf.type="plain")

#PART 1
#Survival probability for American at 10 years
S10A = summary(fitA, times = 10)
Quantity1 = S10A$surv  ; Quantity1
Quantity2 = S10A$lower ; Quantity2
Quantity3 = S10A$upper ; Quantity3

#Probability British to die by 12 years
S12B = summary(fitB, times = 12)
Quantity4 = 1 - S12B$surv  ; Quantity4
Quantity5 = 1 - S12B$upper ; Quantity5
Quantity6 = 1 - S12B$lower ; Quantity6

#Conditional probability American dies between 3 and 8 given alive at 3
S3A = summary(fitA, times=3)$surv
S8A = summary(fitA, times=8)$surv
Quantity7 = (S3A - S8A)/S3A ; Quantity7

#Conditional probability British dies between 3 and 8 given alive at 3
S3B = summary(fitB, times=3)$surv
S8B = summary(fitB, times=8)$surv
Quantity8 = (S3B - S8B)/S3B ; Quantity8


#Population proportions of each breed
piA = mean(Group == 1)
piB = mean(Group == 2)

#Survival at 8.2 years
S82A = summary(fitA, times = 8.2, extend = TRUE)$surv
S82B = summary(fitB, times = 8.2, extend = TRUE)$surv

#Probability of death by 8.2
PA_event = 1 - S82A
PB_event = 1 - S82B

#Bayes theorem, probability dead cat by 8.2 is American
Quantity9 = (piA * PA_event) / (piA * PA_event + piB * PB_event) ; Quantity9


#PART 2
#Indicator variable: 1 = British, 0 = American
z = ifelse(Group==2,1,0)
time = Lifetime
status = Censor

#Interpretation of beta sign (given in question)
Quantity10 = 4 ; Quantity10

#Cox partial log-likelihood function
Quantity11 = function(beta){
  loglik = 0
  for(i in which(status==1)){
    riskset = which(time >= time[i])
    loglik = loglik + z[i]*beta - 
      log(sum(exp(z[riskset]*beta)))
  }
  return(loglik)
} ; Quantity11

#Log-likelihood at beta = 0.1
Quantity12 = Quantity11(0.1) ; Quantity12

#Finding MLE of beta by maximising partial likelihood
opt = optimize(function(b) -Quantity11(b), interval=c(-3,3))
Quantity13 = opt$minimum ; Quantity13

#Likelihood ratio test statistic
LRT = 2*(Quantity11(Quantity13) - Quantity11(0))
Quantity14 = LRT               ; Quantity14
Quantity15 = 1 - pchisq(LRT,1) ; Quantity15


#Decision at 2% significance level
Quantity16 = ifelse(Quantity15 < 0.02, 1, 2) ; Quantity16

#Cox inappropriate since infection changes hazard over time (given in question)
Quantity17 = 3 ; Quantity17


#PART 3
#Log-likelihood for illness-death model
Quantity18 = function(mu, sigma, nu){
  loglik = 0
  
  for(i in 1:length(Lifetime)){
    t = Lifetime[i]
    fip = FIPAge[i]
    d = status[i]
    
    if(fip==0){  # never infected
      if(d==1){
        loglik = loglik + log(mu) - (mu+sigma)*t
      } else{
        loglik = loglik - (mu+sigma)*t
      }
    } else{      # infected at fip
      if(d==1){
        loglik = loglik + log(sigma) - (mu+sigma)*fip +
          log(nu) - nu*(t-fip)
      } else{
        loglik = loglik + log(sigma) - (mu+sigma)*fip -
          nu*(t-fip)
      }
    }
  }
  return(loglik)
} ; Quantity18

#Evaluate log-likelihood at given parameters
Quantity19 = Quantity18(0.1,0.1,0.1) ; Quantity19


infected = which(FIPAge>0)

deaths_inf = sum(status[infected]==1)
time_inf = sum(Lifetime[infected] - FIPAge[infected])

#MLE of nu = deaths while infected / infected exposure time
Quantity20 = deaths_inf / time_inf ; Quantity20

#Hazard derived from survival function (given in question)
Quantity21 = 2 ; Quantity21


a = 0.1
b = 0.08

S = function(x) exp(-a*x - 0.5*b*x^2)

#Conditional probability infected between 7 and 10
Quantity22 = 1 - S(10)/S(7) ; Quantity22

###############################################################################
quantities = c(
  Quantity1, Quantity2, Quantity3, Quantity4, Quantity5, Quantity6,
  Quantity7, Quantity8, Quantity9,
  Quantity10, NA, Quantity12, Quantity13, Quantity14, Quantity15,
  Quantity16, Quantity17,
  NA, Quantity19, Quantity20, Quantity21, Quantity22
)

df_excel = data.frame(
  paste0("Quantity", 1:22),
  quantities
)

write.table(
  df_excel,
  "BASIL_H00435828.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)