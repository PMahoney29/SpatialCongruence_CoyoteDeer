data {
  // Coyote 3rd Order space use data
  int<lower=0> Nc;                        // number of space use records
  //int<lower=1> Yc;                      // number of years (for random effect)
  int<lower=1> Ic;                        // number of individuals (for random effect)
  int<lower=1> Kc;                        // number of predictors (fixed effects)
  int<lower=0, upper=1> y_c[Nc];          // Response
  //int<lower=1,upper=Yc> yearsCoy[Nc];   // The years (must be sequential starting from 1)
  int<lower=1, upper=Ic> indsCoy[Nc];     // The individuals
  matrix[Nc,Kc] xC;                       // Predictors
  
  // Coyote removal data  
  int<lower=0> Nr;                       // number of removal records 
  int<lower=1> Na;                       // number of study areas (for random effect)
  int<lower=1> Kr;                       // number of predictors (fixed effects)
  int<lower=0, upper=1> y_r[Nr];         // Response
  int<lower=1,upper=Na> areas[Nr];       // Study Area (north) effect (must be sequential starting from 1)
  matrix[Nr,Kr] xR;                      // Predictors
  matrix[Nr,Kc] xInt;                    // Predictors for intermediate dataset used in estimating coyote use
  
  // Deer 2nd Order space use data
  int<lower=0> Nd;                        // number of space use records
  int<lower=1> Yd;                        // number of years (for random effect)
  int<lower=1> Kd;                        // number of predictors (fixed effects)
  int<lower=0, upper=1> y_d[Nd];          // Response
  int<lower=1,upper=Yd> yearsDeer[Nd];    // The years (must be sequential starting from 1)
  matrix[Nd,Kd] xD;                       // Predictors
  
  // Data for prediction, and Earth Mover's Distance
  // Estimated externally
  int<lower=1> Npred;                  // Data (size) subset for model predictions and fit using GAM
  matrix[Npred, Kc] xCoyPred;          // Predictors for Deer extrapolation
  matrix[Npred, Kr] xRemPred;          // Predictors for Removal extrapolation
  matrix[Npred, Kd] xDeerPred;         // Predictors for Deer extrapolation
}

parameters {
  // Params for the Coyote Resource Selection
  //real alphaCoy;                // intercept
  //vector[Yc] alphaCoyYear;      // random intercept, year
  vector[Ic] alphaCoyInd;         // random intercept, individuals
  vector[Kc] betaCoy;             // beta coeffs for space use model
  //real<lower=0> sigma_alphaCoyYear;
  //real mu_alphaCoyYear;
  
  // Hyperparams for Coyote Resource Selection
  real<lower=0> sigma_alphaCoyInd;
  real mu_alphaCoyInd;
  
  // Params for the removal model
  //real alphaRem;                 // intercept
  vector[Na] alphaRemArea;         // random intercept, study area
  vector[Kr] betaR;                // beta coeffs for removal model
  real betaCoyUse;                 // beta for Use derived from xUse, intermediate data
  
  // Hyperparams for Coyote Removal Risk
  real<lower=0> sigma_alphaRemArea;
  real mu_alphaRemArea;
  
  // Params for the Deer space use model
  //real alphaDeer;                // intercept
  vector[Yd] alphaDeerYear;        // random intercept, year
  vector[Kd] betaDeer;             // beta coeffs for space use model
  
  // Hyperparams for Deer Site Selection
  real<lower=0> sigma_alphaDeerYear;
  real mu_alphaDeerYear;
}

transformed parameters {
  vector<lower=0, upper=1>[Nr] xCoyUse;

  // Estimate Coyote use for removal data
  for (nr in 1:Nr) 
    xCoyUse[nr] = inv_logit(mu_alphaCoyInd + xInt[nr]*betaCoy);
}
model {
  vector[Nc] y_c_hat;  // Local parameter for log odds
  vector[Nr] y_r_hat;  // Local parameter for log odds
  vector[Nd] y_d_hat;  // Local parameter for log odds

  // NOTE: the normal scale parameters represent standard deviations in Stan

  // Coyote Resource Selection
    // Priors
  //alphaCoyYear ~ normal(mu_alphaCoyYear, sigma_alphaCoyYear); 
  alphaCoyInd ~ normal(mu_alphaCoyInd, sigma_alphaCoyInd);
  betaCoy ~ normal(0, 5);  

    // Hyperpriors
  //sigma_alphaCoyYear ~ cauchy(0, 2.5); 
  //mu_alphaCoyYear ~ normal(0, 10);   
  sigma_alphaCoyInd ~ cauchy(0, 2.5); 
  mu_alphaCoyInd ~ normal(0, 10); 

    // Likelihood
  for (nc in 1:Nc) 
    y_c_hat[nc] = alphaCoyInd[indsCoy[nc]]; 
  y_c_hat = y_c_hat + xC*betaCoy; 
  y_c ~ bernoulli_logit(y_c_hat);
  
  
  // Coyote Removal Risk
    // Priors for space use model
  alphaRemArea ~ normal(mu_alphaRemArea, sigma_alphaRemArea); 
  //alphaRemArea ~ normal(0, 1); 
  betaR ~ normal(0, 5); 
  betaCoyUse ~ normal(0, 5);
  
    // Hyperpriors
  sigma_alphaRemArea ~ cauchy(0, 2.5); 
  mu_alphaRemArea ~ normal(0, 10);   

    // Likelihood
  for (nr in 1:Nr) 
    y_r_hat[nr] = alphaRemArea[areas[nr]]; 
  y_r_hat = y_r_hat + xR*betaR + xCoyUse*betaCoyUse;
  y_r ~ bernoulli_logit(y_r_hat); 
  
  
  // Deer Resource Selection
    // Priors
  alphaDeerYear ~ normal(mu_alphaDeerYear, sigma_alphaDeerYear); 
  betaDeer ~ normal(0, 5);  
  
    // Hyperpriors
  sigma_alphaDeerYear ~ cauchy(0, 2.5); 
  mu_alphaDeerYear ~ normal(0, 10); 

    // Likelihood
  for (nd in 1:Nd) 
    y_d_hat[nd] = alphaDeerYear[yearsDeer[nd]]; 
  y_d_hat = y_d_hat + xD*betaDeer; 
  y_d ~ bernoulli_logit(y_d_hat);
}

generated quantities{
  // NOTE: I do not recommend running both PPcheck and EMD in the same model run.
  // It will likely over-run your memory for large spatial data sets.
  
  // Posterior predictive check variables, comment out if not needed
  //vector[Nc] y_predC;
  //vector[Nr] y_predR;
  //vector[Nd] y_predD;
  
  // Variables to predict output for Earth Mover's Distance and Congruence
  vector[Npred] x_CoyUsePred;
  vector[Npred] x_RemPred;  
  vector[Npred] y_DeerPred;
  
  // Posterior predictive check, comment out if not needed
  //for (nc in 1:Nc) 
  //  y_predC[nc] = bernoulli_rng(inv_logit(alphaCoyInd[indsCoy[nc]] + xC[nc]*betaCoy));
  //
  //for (nr in 1:Nr) 
  //  y_predR[nr] = bernoulli_rng(inv_logit(alphaRemArea[areas[nr]] + xR[nr]*betaR + xCoyUse[nr]*betaCoyUse));
  //  
  //for (nd in 1:Nd) 
  //  y_predD[nd] = bernoulli_rng(inv_logit(alphaDeerYear[yearsDeer[nd]] + xD[nd]*betaDeer));
    
  // Output for Earth Mover's Distance and Congruence
  for (nd in 1:Npred) 
    x_CoyUsePred[nd] = inv_logit(mu_alphaCoyInd + xCoyPred[nd]*betaCoy);
  x_RemPred = mu_alphaRemArea + xRemPred*betaR + x_CoyUsePred*betaCoyUse; 
  y_DeerPred = mu_alphaDeerYear + xDeerPred*betaDeer;
}
