data {
    int<lower=0> nrep;
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> nsp;
    int<lower=0> rep[n];
    int<lower=0> year[n];
    real count[n];
    int<lower=0> sp[n];
}


parameters {
    matrix[nyear,nsp] a;
    real b[nrep];
    real <lower=0.00001> sigb;
    real <lower=0> sig;
}

transformed parameters {
   vector[n] mu;
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]]+b[rep[i]]);
}
}

model {
sigb~gamma(100,100);
sig~gamma(100,100);
sig~gamma(100,100);
to_vector(b)~normal(0,sigb);
to_vector(a)~normal(0,100);

count~lognormal(mu,sig);
}

generated quantities { 
vector[n] newy;
matrix[nyear-1,nsp] r;
for(i in 1:n)
    {
newy[i]=lognormal_rng(mu[i],sig);
    }
    
for(j in 1:(nyear-1))
{
  for(k in 1:nsp)
  {
  r[j,k]=a[j+1,k]-a[j,k];
  }
}
}
