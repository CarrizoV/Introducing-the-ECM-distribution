# R functions for the ECM distribution

The following `R` files contain functions required to reproduce results in Carrizo Vergara et al (2025) and can be used for doing further analysis.

The file `ECM_Gauss_ll.R` contains a function to compute the Gaussian pseudo-likelihood for an ECM count arrangement (list of vectors) with given size and (two-times) path probabilities.

The file `ECM_Pairwise_ll.R` contains a function to compute the pairwise composite likelihood for an ECM count arrangement with given size and (two-times) path probabilities.

The file `ECM_Poisson_Gauss_ll.R` contains a function to compute the Gaussian pseudo-likelihood for an ECM-Poisson count arrangement with given size rate and (two-times) path probabilities.

The file `ECM_Poisson_Pairwise_ll.R` contains a function to compute the pairwise composite likelihood for an ECM-Poisson count arrangement with given size rate and (two-times) path probabilities.

The file `dBivPois.R` contains a vectorized function to compute bivariate Poisson log-likelihoods.

The file `Mean_Cov_Functions.R` contains functions to compute mean vectors and covariance matrices associated to some Gaussian trajectories.

The file `P_Snapshot.R` contains a function which, among other utilities, computes one-time and two-times path probabilities for diverse underlying movement trajectory models.

The file `Sim_Brown.R` contains a function to simulate Brownian trajectories.

The file `Sim_OU.R` contains a function to simulate Ornstein-Uhlenbeck trajectories.

The file `Sim_Population_Count_list.R` contains a function to simulate survey counts of continously moving individuals according to some continuous-time trajectory model over possibly irregularly sampled squares.

The file `Voting_ll_Gauss.R` contains function to compute the Gaussian pseudo-likelihood of transfer proportions for Poisson-multinomial counts modeling two-rounds election results. It also contains a function to compute its gradient and convenient softmax-transformation functions.
