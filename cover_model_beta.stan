data {
    int<lower=0> nrep;
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> nsp;
    int<lower=0> rep[n];
    int<lower=0> year[n];
    real count[n];
    int<lower=0> sp[n];
    int<lower=0, upper=1> is_discrete[n]; 
    int<lower=-1, upper=1> y_discrete[n];
}


parameters {
    matrix[nyear,nsp] a;
    matrix[nyear,nsp] a2;
    real b[nrep];
    real b2[nrep];
    real <lower=0.00001> sigb;
    real<lower=0,upper=1> phi;
}

transformed parameters {
   real mu[n];
   real mu2[n]; #for 0 and 1
for(i in 1:n)
{
mu[i]=inv_logit(a[year[i],sp[i]]+b[rep[i]]);
mu2[i]=inv_logit(a2[year[i],sp[i]]+b2[rep[i]]);
}
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
to_vector(b)~normal(0,sigb);
phi ~ beta(1, 1);

for(i in 1:n)
{
  if (is_discrete[i] == 1) 
    {
      y_discrete[i] ~ bernoulli(mu2[i]);
    } 
    else 
    {
      count[i]~beta(mu[i]*phi,(1-mu[i])*phi);
    }
}
}

generated quantities { 
real newy[n];
real newy2[n];
matrix[nyear-1,nsp] r;
matrix[nyear-1,nsp] lambda;
matrix[nyear-1,nsp] r2;
matrix[nyear-1,nsp] lambda2;
for(i in 1:n)
    {
newy[i]= beta_rng(mu[i]*phi,(1-mu[i])*phi);
newy2[i]= bernoulli_rng(mu2[i]);
    }
    
for(j in 1:(nyear-1))
{
  for(k in 1:nsp)
  {
  r[j,k]=log(inv_logit(a[j+1,k])/inv_logit(a[j,k]));
  lambda[j,k]=inv_logit(a[j+1,k])/inv_logit(a[j,k]);
  r2[j,k]=log(inv_logit(a2[j+1,k])/inv_logit(a2[j,k]));
  lambda2[j,k]=inv_logit(a2[j+1,k])/inv_logit(a2[j,k]);
  }
}
}


