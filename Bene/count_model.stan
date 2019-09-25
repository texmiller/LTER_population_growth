data {
    int<lower=0> nrep;
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> nsp;
    int<lower=0> rep[n];
    int<lower=0> year[n];
    int<lower=0> count[n];
    int<lower=0> sp[n];
}


parameters {
    matrix[nyear,nsp] a;
    real b[nrep];
    real <lower=0.00001> sigb;
    //real <lower=0.00001> epsilon[nsp]; // species specific epsilon?
    real <lower=0.00001> epsilon;
}

transformed parameters {
   real mu[n];
   real phi[n];
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]] + b[rep[i]]);
//phi[i]=mu[i]+(mu[i]^2)/epsilon[sp[i]]; // species specific epsilon?
phi[i]=mu[i]+(mu[i]^2)/epsilon;
}
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
//to_vector(epsilon)~gamma(100,100); //species specific epsilon?
epsilon~gamma(100,100);
to_vector(b)~normal(0,sigb);

for(i in 1:n)
    {
    count[i]~neg_binomial_2(mu[i],phi[i]); //negative binomial
    }
}

generated quantities { 
vector[n] newcount;
matrix[nyear-1,nsp] r;
for(i in 1:n)
    {
newcount[i]=neg_binomial_2_rng(mu[i],phi[i]);
    }

for(j in 1:(nyear-1))
  {
  for(k in 1:nsp)
    {
    r[j,k]=log(exp(a[j+1,k])/exp(a[j,k]));
    }
  }
}


