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
}

transformed parameters {
   real mu[n];
for(i in 1:n)
{
mu[i]=inv_logit(a[year[i],sp[i]]+b[rep[i]]);
}
}

model {

to_vector(a)~normal(0,100);
sigb~gamma(100,100);
to_vector(b)~normal(0,sigb);

count~binomial(100,to_vector(mu));

}

generated quantities { 
int<lower=0> newy[n];
matrix[nyear-1,nsp] r;
for(i in 1:n)
    {
newy[i]= binomial_rng(100,mu[i]);
    }
    
for(j in 1:(nyear-1))
{
  for(k in 1:nsp)
  {
  r[j,k]=log(inv_logit(a[j+1,k]))-log(inv_logit(a[j,k]));
  }
}
}


