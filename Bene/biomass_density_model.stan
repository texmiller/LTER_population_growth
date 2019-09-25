data {
    int<lower=0> nrep ;
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
    real <lower=0.00001> sig;
}

transformed parameters {
   vector[n] mu;
   vector[n] alpha;
   vector[n] beta;
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]]+b[rep[i]]);
}
alpha=mu .*mu/sig;
beta=mu/sig;
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
sig~gamma(100,100);
to_vector(b)~normal(0,sigb);

count~gamma(alpha,beta);
}

generated quantities { 
vector[n] newy;
matrix[nyear-1,nsp] r;
for(i in 1:n)
    {
newy[i]=gamma_rng(alpha[i],beta[i]);
    }
    
for(j in 1:(nyear-1))
{
  for(k in 1:nsp)
  {
  r[j,k]=log(exp(a[j+1,k])/exp(a[j,k]));
  }
}
}


