data {
    int<lower=0> nrep;
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> rep[n];
    int<lower=0> year[n];
    int<lower=0> count[n];
}


parameters {
    real a[nyear];
    real b[nrep];
    real <lower=0.00001> sigb;
}

transformed parameters {
   real mu[n];
for(i in 1:n)
{
mu[i]=exp(a[year[i]]+b[rep[i]]);
}
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
to_vector(b)~normal(0,sigb);

for(i in 1:n)
    {
    count[i]~poisson(mu[i]);
    }
}

generated quantities { 
vector[n] newcount;
for(i in 1:n)
    {
newcount[i]=poisson_rng(mu[i]);
    }
}


