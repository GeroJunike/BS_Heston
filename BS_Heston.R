########### Black-Scholes model ###########
#Characteristic function of the centralized log-returns in the BS model. params is volatility
phi = function(u, mat, params, S0, r){
  return(exp(-params^2 * mat * u^2 / 2))
}

#mu is equal to E[log(S_T)]
mu = function(mat, params, S0, r){
  return(log(S0) + r * mat - 1/2 * params^2 * mat)
}

#cosine coefficients of the density.
ck = function(L, mat, N, params, S0, r){
  k = 0:N
  return(1 / L * Re(phi(k * pi / (2 * L), mat, params, S0, r) * exp(1i * k * pi/2)))
}

#cosine coefficients of a put option, see Appendix Junike and Pankrashkin (2022).
vk = function(K, L, mat, N, params, S0, r){
  mymu = mu(mat, params, S0, r) #mu = E[log(S_T)]
  d = min(log(K) - mymu, L)
  if(d <= -L)
    return(rep(0, N + 1)) #Return zero vector
  k = 0:N
  psi0 = 2 * L / (k * pi) * (sin(k * pi * (d + L) / (2 * L)))
  psi0[1] = d + L
  tmp1 = k * pi / (2 * L) * sin( k * pi * (d + L) / (2 * L))
  tmp2 = cos(k * pi * (d + L) / (2 * L))
  tmp3 = 1 + (k * pi / (2 * L))^2
  psi1 = (exp(d) * (tmp1 + tmp2) - exp(-L)) / tmp3
  return(exp(-r * mat) * (K * psi0 - exp(mymu) * psi1))
}

#approximation of put option by COS method
put_COS = function(K, L, mat, N, params, S0, r){
  tmp = ck(L, mat, N, params, S0, r) * vk(K, L, mat, N, params, S0, r)
  tmp[1] = 0.5 * tmp[1] #First term is weighted by 1/2
  return(sum(tmp))
}

#approximation of call option by COS method using put-call parity
call_COS = function(K, L, mat, N, params, S0, r){
  return(put_COS(K, L, mat, N, params, S0, r) + S0 - K * exp(-r * mat))
}

#Price a call option in BS model by the COS method
eps = 10^-4 #error tolerance
K = 90 #strike
S0 = 100 #current stock price
r = 0.1 #interest rates
params = 0.2 #volatility
mat = 0.7 #maturity
n = 8 #number of moments to get truncation range
mu_n = (params * sqrt(mat))^n * prod(seq(1, n, by = 2)) #n-th moment of log-returns
L = (2 * K * exp(-r * mat) * mu_n / eps)^(1 / n) #Truncation range, Junike (2024, Eq. (3.10))
s = 40 #number of derivatives to determine the number of terms
boundDeriv = gamma(s / 2 + 1) * 2^(s / 2) / (pi * (params * sqrt(mat))^(s + 2))
tmp = 2^(s + 5 / 2) * boundDeriv * L^(s + 2) * 12 * K * exp(-r * mat)
N = ceiling((tmp / (s * pi^(s + 1) * eps))^(1 / s)) #Junike (2024, Sec. 6.1)
call_COS(K, L, mat, N, params, S0, r) #The price of call option is 17.24655.


########### Heston model ###########
#Use functions ck, vk, put_COS, call_COS from the Black-Scholes model above.

#Characteristic function of log-returns in the Heston with parameters params.
#The characteristic function is taken from Schoutens et. al (2004).
psiLogST_Heston = function(u, mat, params, S0, r){
  kappa = params[1] #speed of mean reversion
  theta = params[2] #level of mean reversion
  xi = params[3] #vol of vol
  rho = params[4] #correlation vol stock
  v0 = params[5] #initial vol
  d = sqrt((rho * xi * u * 1i - kappa)^2 - xi^2 * (-1i * u - u^2))
  mytmp = kappa - rho * xi * u * 1i
  g = (mytmp - d) / (mytmp + d)
  expdmat = exp(-d * mat)
  tmp0 = 1i * u * (log(S0) + r * mat)
  tmp1 = (mytmp - d) * mat - 2 * log((1 - g * expdmat) / (1 - g))
  tmp2 = theta * kappa * xi^(-2) * tmp1
  tmp3 = v0 * xi^(-2) * (mytmp - d) * (1 - expdmat) / (1 - g * expdmat)
  exp(tmp0 + tmp2 + tmp3)
}

#install.packages("Deriv")
library(Deriv) #There are much faster alternatives like SageMath.
psiLogST_Heston1=Deriv(psiLogST_Heston, "u")

#mu is equal to E[log(S_T)]
mu = function(mat, params, S0, r){
  Re(-1i * psiLogST_Heston1(0, mat, params, S0, r))
}

#Characteristic function of centralized log-returns in the Heston model.
phi = function(u, mat, params, S0, r){
  psiLogST_Heston(u, mat, params, S0, r) * exp(-1i * u * mu(mat, params, S0, r))
}

#Derivatives of the characteristic function of the centralized log-returns in the Heston model.
if(F){
  phi1 = Deriv(phi, "u")
  phi2 = Deriv(phi1, "u")
  phi3 = Deriv(phi2, "u") #Takes very long but has to be done only once.
  phi4 = Deriv(phi3, "u") #Takes very long but has to be done only once.
  save(phi4, file = "phi4.RData") #save for later use. Load with load("phi4.RData").
}
source("phi4.R") #loads the forth derivative

#Price a put option in Heston model by the COS method.
eps = 10^-1 #error tolerance
K = 90 #strike
S0 = 100 #current stock price
r = 0.1 #interest rates
params = c(0.6067, 0.0707, 0.2928, -0.7571, 0.0654)
mat = 0.7 #maturity
mu_n = abs(phi4(0, mat, params, S0, r)) #4-th moment of log-returns.
L = (2 * K * exp(-r * mat) * mu_n / eps)^(1 / 4) #Junike (2024, Eq. (3.10)).
s = 19 #number of derivatives to determine the number of terms
integrand = function(u){1 / (2 * pi) * abs(u)^(s + 1) * abs(phi(u, mat, params, S0, r))}
boundDeriv = integrate(integrand, -Inf, Inf)$value
tmp = 2^(s + 5 / 2) * boundDeriv * L^(s + 2) * 12 * K * exp(-r * mat)
N = ceiling((tmp / (s * pi^(s + 1) * eps))^(1 / s)) #Number of terms, Junike (2024, Sec. 6.1)
put_COS(K, L, mat, N, params, S0, r) #The price of put option is 2.773954.
