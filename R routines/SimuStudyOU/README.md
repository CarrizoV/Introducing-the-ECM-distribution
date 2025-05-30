# Instructions


The `R` routines `SimuStudyOUstatToy2D(parallel).R` and `SimuStudyOUstatToy2DUnknownN(parallel).R` contain code to launch simulation studies. Both simulate underlying steady-state Ornstein-Uhlenbeck individuals, count them over survey squares and time defined in `ToySetting2D.R`, and then use MGLE and/or MCLE method for estimation. The user can propose values for the activity center $(z_{1},z_{2})$, the home range $\tau$ and the infinitesimal quadratic speed $\sigma$. In `SimuStudyOUstatToy2D(parallel).R` the number of individuals is known (ECM distribution), in `SimuStudyOUstatToy2DUnknownN(parallel).R` the number of individuals is unknown and assumed Poisson (ECM-Poisson distribution). The parameter setting used in the paper $(\tau,\sigma,z_{1},z_{2}) = ( 0.4 , 0.4\sqrt{0.002} , -0.2 , 0.1 )$  is set as default value.  The codes are parallelized, `set.seed` has been applied so the same results as in the paper can be reproduced, although remark the need of a server with `95` cores for it. For each parameter setting, the computing time of `1045` simulation-estimations obtained by doing `11` simulation-estimations in each of the `95` cores took around 6 hours. As noted in Appendix B in Carrizo Vergara et al. (2025), the case of the MGLE for ECM-Poisson with $\lambda = 10^2$ requires extra changes. Those are explained in the same code, in the comments.


The results obtained by us are already present in the folder `output`. The routine `Analysis.R` is used to read them. This routine will search in the output folder for the results of a simulation study with a given parameter settings for different quantity of individuals and/or size rates at once. It will take out erratic estimations (see Appendix B in paper) retaining `1000` samples for analysis. It will then produce comparative boxplots and correlograms for the samples. It will also compute diverse bias statistics.


## References

Carrizo Vergara, R., Kéry, M., & Hefley, T. (2025). The evolving categories multinomial distribution: introduction with applications to movement ecology and vote transfer. arXiv preprint 2505.20151.


