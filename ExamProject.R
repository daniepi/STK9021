#### STK 9021 - Project Exam
## Exercise 1
## Cooling of new-born babies. Two groups {non-cooled, cooled}
## y0 ~ Binom(m0, theta0); theta0_MLE = 22/79
## y1 ~ Binom(m1, theta1); theta1_MLE = 19/78
## Of interest rho=theta1/theta0 - risk ratio
## a.) Jeffreys for binomial is Beta(1/2,1/2)
## b.) This corresponds to Beta posterior Beta(1/2+y, 1/2+m-y)
library("ggplot2", lib.loc="~/R/win-library/3.4")
library("tikzDevice", lib.loc="~/R/win-library/3.4")
setwd("D:/School/PhD/Fall2017/STK9021/Project")

babies = data.frame(dd = c(22,19), m=c(79,78))
theta_ML = babies$dd/babies$m
n.sim=10**5; alpha=beta=0.5
theta0_post = rbeta(n.sim, alpha+babies$dd[1], beta+babies$m[1]-babies$dd[1])
theta1_post = rbeta(n.sim, alpha+babies$dd[2], beta+babies$m[2]-babies$dd[2])
rho_post = theta1_post/theta0_post
#hist(rho_post, freq=F)
toplot = data.frame(X=rho_post)
#tikz('JeffreysPost.tex', height=4.5, width=6.5)
ggplot(toplot, aes(X)) +
  stat_bin(aes(y=..count../sum(..count..)), col='black', fill='white') +
  #geom_histogram(col='black', fill='white')+
  labs(y = "Density", x = 'rho')+
  theme(panel.background = element_rect(fill=alpha('lightblue', 0.3)))
#dev.off()
round(quantile(rho_post, probs = c(0.025,0.5,0.975)),4)
##
theta_seq = seq(0,1,0.01)
toplot=data.frame(y=dbeta(theta_seq, 0.5, 0.5), x=theta_seq)
#tikz('JeffreysDens.tex', height=4.5, width=6.5)
ggplot(data = toplot, aes(x = x, y = y)) +
  geom_line(lwd=1.2) + 
  labs(x = "theta", y = 'jeff')+
  theme(panel.background = element_rect(fill=alpha('lightblue', 0.3)))
#dev.off()
#plot(theta_seq, dbeta(theta_seq, 0.5, 0.5),'l', xlab='theta',ylab='pjeff')
## d.) MLE
rho_hat = theta_ML[2]/theta_ML[1]
gamma_hat = log(rho_hat)
SE_gamma = sqrt((1/babies$m[2])*(1/theta_ML[2] -1)+(1/babies$m[1])*(1/theta_ML[1] -1))
## 95% CI
gamma_hat + qnorm(c(0.025,0.975))*SE_gamma ## gamma
exp(gamma_hat + qnorm(c(0.025,0.975))*SE_gamma) ## rho
### g.) Normal prior for gamma ~ N(0, 0.35**2) ==> rho ~ LN(0, 0.35**2)
## Find prior 95% CI
sigma_0 = 0.35; gamma_0 = 0
gamma_0 + qnorm(c(0.025,0.975))*sigma_0
exp(gamma_0 + qnorm(c(0.025,0.975))*sigma_0)
rho_LN = exp(gamma_0+sigma_0**2/2)
var_LN = rho_LN**2*(exp(sigma_0**2)-1)
rho_sim = rlnorm(n.sim*2, gamma_0, sigma_0)
quantile(rho_sim, probs = c(0.025,0.5,0.975))
## h.) Plot to one figure Log-normal prior, posterior from Beta priors (Jeffreys) and posterior
## from log-normal prior.
gamma_B = gamma_0 + sigma_0**2/(sigma_0**2 + SE_gamma**2)*(gamma_hat-gamma_0)
sigma2_B  = (sigma_0**2*SE_gamma**2)/(sigma_0**2+SE_gamma**2)
gamma_sim = rnorm(n.sim, gamma_B, sqrt(sigma2_B))
rho_sim = rlnorm(n.sim, gamma_B, sqrt(sigma2_B))
### Plot
rho_val = seq(0,4,0.01)
plot(rho_val, dlnorm(rho_val, 0, sigma_0),'l', xlab='rho',ylab='f', ylim=c(0,2.05))
lines(rho_val, dlnorm(rho_val, gamma_B, sqrt(sigma2_B)),col='red')
lines(density(rho_post, n=length(rho_val)), col='blue')

toplot = data.frame(y= c(dlnorm(rho_val, 0, sigma_0), dlnorm(rho_val, gamma_B, sqrt(sigma2_B)),
                         density(rho_post, n=length(rho_val))$y), x=c(rep(rho_val,2), density(rho_post, n=length(rho_val))$x), 
          Est=rep(c('LN prior', "LN posterior", "Jeffreys posterior"), each=length(rho_val)))
#tikz('rho_post.tex', height=4.5, width=6.5)
ggplot(data = toplot, aes(x = x, y = y, col=Est)) +
  geom_line(aes(linetype=Est, size=Est)) + 
  scale_linetype_manual(values=c("dashed", "dashed",'solid'))+
  scale_size_manual(values=c(1,1, 1.2))+
  labs(x = "rho", y = 'f') +
  theme(panel.background = element_rect(fill=alpha('lightblue', 0.3)), legend.title = element_blank())
#dev.off()
### Quantiles
round(mean(rho_post),4)
round(quantile(rho_post, probs=c(0.025,0.5,0.975)),4)
round(mean(rho_sim),4); round(exp(gamma_B + sigma2_B/2),4)
round(quantile(rho_sim, probs=c(0.025,0.5,0.975)),4)

mean(rho_post<1); mean(rho_sim<1)
## i.) P(theta1 < theta0 | data) = ?
mean(theta1_post < theta0_post)
## For what range of possibe outcomes y1 in {0,m1} is this probability 0.95 or higher?
theta0_post = rbeta(n.sim, alpha+babies$dd[1], beta+babies$m[1]-babies$dd[1])
theta1_post = rbeta(n.sim, alpha+babies$dd[2]-10, beta+babies$m[2]-babies$dd[2]-10)
mean(theta1_post < theta0_post)
## j.) find Beta priors for theta0 theta1 such that gamma=log(theta1)-log(theta0) has 0 mean and
## 0.35**2 variance.
optimize_beta =  function(params){
  alpha0 = params[1]; beta0 = params[2];
  alpha1 = params[3]; beta1 = params[4]
  (digamma(alpha1)-digamma(alpha1+beta1)+digamma(alpha0+beta0)-digamma(alpha0))**2 +
    (trigamma(alpha1)-trigamma(alpha1+beta1)+trigamma(alpha0)-trigamma(alpha0+beta0)-0.35**2)**2
}
opt = optim(fn=optimize_beta, par=c(3,3,3,3))
alpha0_opt = opt$par[1]; beta0_opt = opt$par[2]
alpha1_opt = opt$par[3]; beta1_opt = opt$par[4]

theta0_sim = rbeta(n.sim*5, alpha0_opt, beta0_opt)
theta1_sim = rbeta(n.sim*5, alpha1_opt, beta1_opt)
gamma=log(theta1_sim/theta0_sim)
mean(gamma); sd(gamma)

theta0_post = rbeta(n.sim, alpha0_opt+babies$dd[1], beta0_opt+babies$m[1]-babies$dd[1])
theta1_post = rbeta(n.sim, alpha1_opt+babies$dd[2], beta1_opt+babies$m[2]-babies$dd[2])
rho_post = theta1_post/theta0_post
hist(rho_post)


### Exercise 2
### Dosage level
### 10 levels - in each level Binomial experiment with m=5, y_j ~ Binom(m, theta_j)
### Where we model theta_j ~ logit^-1(a+bx_j) = exp(a+bx_j)/(1+exp(a+bx_j))
## x_j ~ dosage levels
rats = scan('clipboard', sep=',')
x=scan('clipboard', sep=',')
m=5
theta_MLE = rats/m
## Binomial log-likelihood
logistic = function(x, params){
  a=params[1]; b=params[2]
  1/(1+exp(-a-b*x))
}
binom_loglik = function(y,x,params){
  a=params[1]; b=params[2]
  theta = logistic(x, params)
  sum(y*log(theta)+(m-y)*log(1-theta))
}
opt = optim(binom_loglik, par=c(1,1), y=rats, x= x, hessian = T, control=list(fnscale=-1))
a_MLE = opt$par[1]; b_MLE = opt$par[2]
Hessian = opt$hessian
Sigma = solve(-Hessian)
rho=Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2])
SE_MLE = sqrt(diag(Sigma))
Chol = chol(Sigma)
eps = matrix(rnorm(2*n.sim),ncol=2)
ab_Normal = t(c(a_MLE, b_MLE) + t(Chol)%*%t(eps))
### b.) Simulate from posterior for (a,b) via MCMC, where a,b are flat on [-10,10]
n.sim=10**4
MH_rats = function(n.sim, prior.a = c(-10,10), prior.b=c(-10,10)){
  a=b=rep(0, n.sim)
  prob1 = runif(n.sim-1); prob2=runif(n.sim-1)
  acc.a=acc.b=0
  for(i in 2:n.sim){
    a.prop = rnorm(1, a[i-1], 1); b.prop = rnorm(1, b[i-1], 1)
    if(a.prop > prior.a[1] & a.prop < prior.a[2]){
      lnew = binom_loglik(rats, x, params=c(a.prop, b[i-1]))
      lold = binom_loglik(rats, x, params=c(a[i-1], b[i-1]))
      if(prob1[i-1] < exp(lnew-lold)){
        a[i] = a.prop
        acc.a = acc.a+1
      } else{
        a[i] = a[i-1]
      }
    } else{
      a[i] = a[i-1]
    }
    if(b.prop > prior.b[1] & b.prop < prior.b[2]){
      lnew = binom_loglik(rats, x, params=c(a[i], b.prop))
      lold = binom_loglik(rats, x, params=c(a[i], b[i-1]))
      if(prob2[i-1] < exp(lnew-lold)){
        b[i] = b.prop
        acc.b = acc.b+1
      } else{
        b[i] = b[i-1]
      }
    } else{
      b[i] = b[i-1]
    }
  }
  burnin = n.sim/10
  a=a[(burnin+1):n.sim]
  b=b[(burnin+1):n.sim]
  return(list(a=a,b=b,acc.a=acc.a/(n.sim-1),acc.b=acc.b/(n.sim-1)))
}
tmp=MH_rats(n.sim)
a=tmp$a; b=tmp$b
plot.ts(a); plot.ts(b)
acf(a); acf(b)
a_MCMC =mean(a); SE_a_MCMC = sd(a)
b_MCMC =mean(b); SE_b_MCMC = sd(b)
rho_MCMC = cor(a,b)
x_val = seq(0,2.5,0.01)
fit=mapply(FUN=function(a,b){logistic(x_val, c(a,b))}, a, b)
fitted = apply(fit, 1, mean)
q_low = apply(fit, 1, quantile, probs=0.05)
q_high = apply(fit, 1, quantile, probs=0.95)
### Normal approximation
fit=mapply(FUN=function(a,b){logistic(x_val, c(a,b))}, ab_Normal[,1], ab_Normal[,2])
q_lowN = apply(fit, 1, quantile, probs=0.05)
q_highN = apply(fit, 1, quantile, probs=0.95)
## Plot
#tikz('fitted.tex', height=4.5, width=8.5)
plot(x, theta_MLE, xlim=c(0,2.5), ylim=c(-0.1,1), pch=16)
lines(x_val, logistic(x_val, c(a_MLE, b_MLE)), lwd=2)
lines(x_val, fitted, lwd=2, col='red', lty=3)
lines(x_val, q_low, lwd=2, col='blue',lty=2)
lines(x_val, q_high, lwd=2, col='blue',lty=2)
#lines(x_val, q_lowN, lwd=2, col='green',lty=2)
#lines(x_val, q_highN, lwd=2, col='green',lty=2)
legend('topleft',bty='n',lty=c(1,3,2,2),col=c(1,2,'blue','blue'),lwd=rep(2,4),
       legend=c('thetaMLE', 'thetaMCMC','5','95'),ncol=2)
legend('bottomright',bty='n',pch=16,legend='Raw')
#dev.off()
## d.) Posterior histogram for theta_j = P(A| x=2.5), i.e. new dosage level
x_new=mapply(FUN=function(a,b){logistic(2.5, c(a,b))}, a, b)
hist(x_new, freq=F); lines(density(x_new),col='red')
toplot = data.frame(X=x_new)
#tikz('posterior_new.tex', height=4.5, width=6.5)
ggplot(toplot, aes(X)) +
  stat_bin(aes(y=..density..), col='black', fill='white') +
  geom_density(alpha=.2, fill='#FF6666')+
  #geom_histogram(col='black', fill='white')+
  labs(y = "Density", x = 'theta')+
  theme(panel.background = element_rect(fill=alpha('lightblue', 0.3)))
#dev.off()
## e.) Predictive for y_new ~ Binom(m=5, theta_new)
y_new = rbinom(length(x_new),size=m, prob=x_new)
y_new_med = qbinom(0.5, m, prob=x_new)
toplot=data.frame(x=rep(0:5, times=2), y=c(prop.table(table(y_new)), prop.table(table(y_new_med))), fill=c(rep(c('ynew','ymed'), each=6)))

#tikz('posterior_y.tex', height=4.5, width=6.5)
ggplot(toplot, aes(x=x,y=y, fill=fill)) +
  geom_col(position='dodge', width=0.5) +
  labs(x = "y", y = 'Density') +
  theme(panel.background = element_rect(fill=alpha('lightblue', 0.3)), legend.title = element_blank())+
  scale_x_continuous(breaks=0:5,labels=0:5)
#dev.off()


### f.) new prior p(a,b) ~ flat on [-10,10]x[0,10] i.e. only upward trend.
tmp = MH_rats(n.sim, prior.b = c(0,10))
tmp$acc.a; tmp$acc.b
a=tmp$a; b=tmp$b
plot.ts(a); plot.ts(b)
acf(a); acf(b)
### Posterior for LD-50 parameter
x0 = (log(1)-a)/b
round(quantile(x0, probs=c(0.05,0.5,0.95)),4)
