data {
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> nsp;
    int<lower=0> year[n];
    real count[n];
    int<lower=0> sp[n];
}


parameters {
    matrix[nyear,nsp] a;
    real<lower=0.001,upper=100> sig;
}

transformed parameters {
   vector[n] mu;
   vector[n] alpha;
   vector[n] beta;
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]]);
}
alpha=square(mu)/sig;
beta=mu/sig;
}

model {

to_vector(a)~cauchy(0,2.5);
sig~uniform(0,100);

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
  r[j,k]=a[j+1,k]-a[j,k];
  }
}
}


