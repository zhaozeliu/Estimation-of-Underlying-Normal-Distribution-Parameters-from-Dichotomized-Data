
library(bin2norm)  # Custom package
library(optparse)  # Command-line argument parser
library(lme4)      # Mixed models
library(rstan)     # Bayesian inference
library(statmod)
library(readr)
ukb35193 <- read.csv("C:/Users/78663/liu/lzz/2024-12-14/likelihood/ukb35193.csv")
# dim 502527 X 18 first column sample id 502527 second column 22 total study numbers

realdataanalysis=function(ukb35193,colnumber)
{

dataukb35193=na.omit(ukb35193[,c(1,2,colnumber)])
varname=colnames(dataukb35193)[3]
par(mfrow = c(1, 2))
qqnorm(dataukb35193[,3],main = paste("Real data QQ Plot"))
qqline(dataukb35193[,3], col = "red", lwd = 2)

hist(dataukb35193[,3], breaks = 30, probability = TRUE, col = "lightblue", main = paste("Real data histogram"),xlab = "data value")
lines(density(dataukb35193[,3]), col = "blue", lwd = 2)  # density curve
curve(dnorm(x, mean = mean(dataukb35193[,3]), sd = sd(dataukb35193[,3])), col = "red", lwd = 2, add = TRUE)  # normal curve


par(mfrow = c(1, 1))

library(moments)
skew=skewness(dataukb35193[,3])
kur=kurtosis(dataukb35193[,3])

dataeachstudysmean=c()
dataeachstudyssd=c()
for (j in 1:length(table(dataukb35193[,2]))) {
  studyjmean=mean(dataukb35193[dataukb35193$X54.0.0==names(table(dataukb35193$X54.0.0))[j],3])
  studyjsd=sd(dataukb35193[dataukb35193$X54.0.0==names(table(dataukb35193$X54.0.0))[j],3])
  dataeachstudysmean=c(dataeachstudysmean,studyjmean)
  dataeachstudyssd=c(dataeachstudyssd,studyjsd)
  
}

###fix mean
i=c()
n_i=c()
c_i=c()
p_i_obs=c()
set.seed(2)
for (j in 1:length(table(dataukb35193[,2]))) {
  i=c(i,j)
  #the j th study data
  studydata=dataukb35193[dataukb35193$X54.0.0==names(table(dataukb35193$X54.0.0))[j],3]
  n_i=c(n_i,length(studydata))
  #Randomly pick one element from the vector as cut point cj
  cutpoint <- sample(studydata, 1)
  c_i=c(c_i,cutpoint)
  #Count how many elements are greater than this value   kj
  count_greater <- sum(studydata > cutpoint)
  #p_j_obs
  frequence=count_greater/length(studydata)
  p_i_obs=c(p_i_obs,frequence)
}

fittedFIXmeandata=data.frame(i,n_i,c_i,p_i_obs)
colnames(fittedFIXmeandata)=c("i","n_i","c_i","p_i_obs")
fixmlemu=estimate_singleThresh_MLE(n_i, c_i, p_i_obs, use_wols_init = TRUE)$mu
fixmlesigma=estimate_singleThresh_MLE(n_i, c_i, p_i_obs, use_wols_init = TRUE)$sigma
fixprobitmu=estimate_singleThresh_probit(n_i, c_i, p_i_obs)$mu
fixprobitsigma=estimate_singleThresh_probit(n_i, c_i, p_i_obs)$sigma

par(mfrow = c(1, 2))

library(ggplot2)
x <- dataeachstudysmean
df <- data.frame(value = dataeachstudysmean)

fixmeanmu=ggplot(df, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue", width = 0.3) +
  geom_hline(yintercept = fixmlemu, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 1.05, y = fixmlemu, label = "fix mle mu", color = "red",
           vjust = -1, hjust = 0, size = 4) +
  geom_hline(yintercept = fixprobitmu, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 0.65, y = fixprobitmu, label = "fix glm mu", color = "red",
           vjust = -1, hjust = 0, size = 4) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue") +
  annotate("text", x = 1.2, y = mean(x),
           label = sprintf("Mean = %.2f\nSD = %.2f", mean(x), sd(x)),
           hjust = 0, vjust = -0.5, size = 3.5) +
  labs(title = "Boxplot for all studys' mean value",
       y = "Value", x = "") +
  theme_minimal()
print(fixmeanmu)


x <- dataeachstudyssd
df <- data.frame(value = dataeachstudyssd)

fixmeansigma=ggplot(df, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue", width = 0.3) +
  geom_hline(yintercept = fixmlesigma, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 1.05, y = fixmlesigma, label = "fix mle sigma", color = "red",
           vjust = -1, hjust = 0, size = 4) +
  geom_hline(yintercept = fixprobitsigma, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 0.65, y = fixprobitsigma, label = "fix glm sigma", color = "red",
           vjust = -1, hjust = 0, size = 4) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue") +
  annotate("text", x = 1.2, y = mean(x),
           label = sprintf("Mean = %.2f\nSD = %.2f", mean(x), sd(x)),
           hjust = 0, vjust = -0.5, size = 3.5) +
  labs(title = "Boxplot for all studys' sd value",
       y = "Value", x = "") +
  theme_minimal()

print(fixmeansigma)
par(mfrow = c(1, 1))
###random mean
data_list <- list(
  n_i = n_i,
  c_ij = vector("list", length(table(dataukb35193[,2]))),
  p_ij_obs = vector("list", length(table(dataukb35193[,2])))
)
set.seed(3)
for(i in 1:length(table(dataukb35193[,2]))){
  studydata=dataukb35193[dataukb35193$X54.0.0==names(table(dataukb35193$X54.0.0))[i],3]
  m_i <- sample(50:100, size=1) #random cutpoints number for each study
  c_ij_vec <- numeric(m_i)
  p_ij_obs_vec <- numeric(m_i)
  
  for(j in seq_len(m_i)){
    c_ij <- sample(studydata,1)
    k_ij <- sum(studydata > c_ij)
    p_ij=k_ij/length(studydata)
    p_ij_obs_vec[j] <- k_ij / n_i[i]
    c_ij_vec[j]<- c_ij
  }
  
  data_list$c_ij[[i]]     <- c_ij_vec
  data_list$p_ij_obs[[i]] <- p_ij_obs_vec
}
fittedRANmeandata=data_list
ranmle=estimate_multiThresh_MLE(fittedRANmeandata)
ranglmm=estimate_multiThresh_GLMM(fittedRANmeandata)
ranmcmc=estimate_multiThresh_MCMC(fittedRANmeandata)

list(fixmlemu=fixmlemu,fixmlesigma=fixmlesigma,fixprobitmu=fixprobitmu,fixprobitsigma=fixprobitsigma,ranmle=ranmle,ranglmm=ranglmm,ranmcmc=ranmcmc,skewness = skew, kurtosis = kur, studymeans=dataeachstudysmean,studysds=dataeachstudyssd)
}


2.661939  #3

0.23^2/(2.62^2+0.23^2)   r0.007  r0.9 7.86
TRYHUGESIZEsmall=generate_multiple_thresholds_data(100, 2.62, 0.007, "A") 
TRYHUGESIZElarge=generate_multiple_thresholds_data(100, 2.62, 0.007, "C") 
estimate_multiThresh_MCMC(TRYHUGESIZElarge)
TRYHUGESIZElarge$n_i=TRYHUGESIZEsmall$n_i

tau <- sqrt((r/(1-r)) * sigma^2)

TRYHUGESIZElarge=generate_multiple_thresholds_data(100, 2.62, 0.007, "A") 

#when A to C, mu0 increasing, tau bad
