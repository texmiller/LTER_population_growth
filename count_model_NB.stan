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
    real phi;
}

transformed parameters {
   real mu[n];
for(i in 1:n)
{
mu[i]=exp(a[year[i],sp[i]] + b[rep[i]]);
}
}

model {
  
phi~cauchy(0,2);
to_vector(a)~normal(0,100);
sigb~gamma(100,100);
to_vector(b)~normal(0,sigb);

for(i in 1:n)
    {
    count[i]~neg_binomial_2(mu[i],phi); //negative binomial
    }
}

//generated quantities { 
//vector[n] newy;
//matrix[nyear-1,nsp] r;
//matrix[nyear-1,nsp] lambda;
//for(i in 1:n)
  // {
//newy[i]=neg_binomial_2_rng(mu[i],phi);
  //  }

//for(j in 1:(nyear-1))
  //{
  //for(k in 1:nsp)
    //{
    //r[j,k]=log(exp(a[j+1,k])/exp(a[j,k]));
    //lambda[j,k]=exp(a[j+1,k])/exp(a[j,k]);
    //}
//  }
//}


