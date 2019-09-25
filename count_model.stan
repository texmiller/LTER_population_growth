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
    real <lower=0.00001> epsilon;
}

transformed parameters {
   real mu[n];
   real phi[n];
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]] + b[rep[i]]);
phi[i]=mu[i]+(mu[i]^2)/epsilon;
}
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
epsilon~gamma(100,100);
to_vector(b)~normal(0,sigb);

for(i in 1:n)
    {
    count[i]~neg_binomial_2(mu[i],phi[i]); #negatice binomial
    }
}

generated quantities { 
vector[n] newcount;
for(i in 1:n)
    {
newcount[i]=neg_binomial_2_rng(mu[i],phi[i]);
    }
}


