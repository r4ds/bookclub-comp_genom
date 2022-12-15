## Chapter 3 CompGenomR
library(matrixStats)
library(mosaic)

#1-2 dataset
set.seed(100)
gset=rnorm(600,mean=200,sd=70)
data=matrix(gset,ncol=6)


#row means/variance
means = rowMeans(data)
vars = rowVars(data)
plot(means)
plot(vars)
boxplot(means)
boxplot(vars)

#SD compared to SE from CLT 
samples = sd(data)
rowSds(data)
clt.se = 70/sqrt(6)  ##standad error = SD/sqrt(n)

data1=matrix(rnorm(6000,mean=200,sd=70),ncol=6)
sd(rowMeans(data1))

data2=matrix(rnorm(6000,mean=200,sd=70),ncol=10)
sd(rowMeans(data2))

#poisson distribution
set.seed(100)
pois1 = rpois(30,5)

quantile((do(1000)*mean(rpois(30,lambda=5)))[,1],probs=c(0.025,0.975))

t.test(pois1)

#### 3.4.3 Linear Models
# set random number seed, so that the random numbers from the text is the same when you run the code.
set.seed(32)

# get 50 X values between 1 and 100
x = runif(50,1,100)

# set b0,b1 and variance (sigma)
b0 = 10
b1 = 2
sigma = 20

# simulate error terms from normal distribution
eps = rnorm(50,0,sigma)

# get y values from the linear equation and addition of error terms
y = b0 + b1*x+ eps

lm(y ~ x)
plot(x,y)
abline(lm(y~x))

summary(lm(y~x))$coefficients[2,4]

###

hmodFile=system.file("extdata",
                     "HistoneModeVSgeneExp.rds",
                     package="compGenomRData")

hdat = readRDS(hmodFile)
plot(hdat$H3k4me3, hdat$measured_log2)
abline(lm(hdat$H3k4me ~ hdat$measured_log2))
summary(lm(hdat$H3k4me ~ hdat$measured_log2))

plot(hdat$H3k27me3, hdat$measured_log2)
abline(lm(hdat$H3k27me ~ hdat$measured_log2))
summary(lm(hdat$H3k27me ~ hdat$measured_log2))

summary(lm(formula = hdat$measured_log2 ~ hdat$H3k4me3 + hdat$H3k27me3))

plot(hdat$H3k27me3, hdat$H3k4me3)
cor(hdat$H3k27me3,hdat$H3k4me3)
