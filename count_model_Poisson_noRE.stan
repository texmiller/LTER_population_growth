data {
    int<lower=0> nyear;
    int<lower=0> n;
    int<lower=0> nsp;
    int<lower=0> year[n];
    int<lower=0> count[n];
    int<lower=0> sp[n];
}


parameters {
    matrix[nyear,nsp] a;
    real <lower=0.00001> phi_od[n];
    real <lower=0.00001> phi;
}

transformed parameters {
   real mu[n];
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]] +phi_od[i]);
}
}

model {

to_vector(a)~normal(0,100);
phi~gamma(100,100); 
to_vector(phi_od)~normal(0,phi);

count~poisson(mu);

}

generated quantities { 
vector[n] newy;
matrix[nyear-1,nsp] r;
for(i in 1:n)
    {
newy[i]=poisson_rng(mu[i]);
    }
for(j in 1:(nyear-1))
{
  for(k in 1:nsp)
  {
  r[j,k]=a[j+1,k]-a[j,k];
  }
}
}


