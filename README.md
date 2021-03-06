# Analysis of the literature_bias dataset
Here you will find our Bayesian analysis on the number of publications (tropical/temparate) in the web of science (WOS) for nine major topics in ecology and evolution.

## For the methods sections:

We performed a generalised linear mixed model with Poisson distribution on the number of papers (N.papers) published in WOS, as explanatory variables, we used REGION (Tropical vs Temparate) and its interactions with TOPIC (see bellow), and TAXON as random effects in Stan (Carpenter 2017) using the R package Rethingking (Mcelreath 2016). To estimate the parameters,  We used Hamilton Monte Carlo (HMC) chaings (4) and 5000 interaction. All non-adaptive priors used were only weakly informative. HMC chains were verified to be well-mixed and stationary. For model comparasion, we use the Widely Applicable Information Criterion (WAIC). WAIC is  a generilized Bayesian version of AIC and has simmilar interpredation  (Watanabe 2010).

## For the results

We can use the figure in the *plot folder* to show the bias in publication numbers given the region and topic, additionally you can find a table with the parameters estimated, and the model comparasion in the *tables folder*  

### Figure 1b (plots/N_papers region by topic.pdf)

Posterior probability density estimates of the number of publications of tropical and temperate publications in major topics of ecology and evolution in the web of science (WOS). The probability of temperate bias in each topic is given at the top-right corner of each panel. Colored shades represent 95% credible intervals (see Table-S# for more details).


### Figure 2 (plots/Ciations region by topic.pdf)
Posterior probability density estimates for the citation rates of tropical and temperate publications in major topics of ecology and evolution in the web of science (WOS). The probability of temperate bias in each topic is given at the top-right corner of each panel. Colored shades represent 95% credible intervals (see Table-S# for more details). 


## Table-S1 (tables/Supplementary Tables.xlsx: plublications) legend

Parameter coefficients estimate from the effects of the region (temperate vs. tropic) and its interactions with sub-topics in ecology and evolution on the number of publications in the web of science. Sub-topics in ecology and evolution are adaptation (AD), climate tolerance (CT), density-dependence (DD), interspecific competition (IC), mimicry (MI), parental care (PC), sexual selection (SS), and speciation (SP). Values of the R statistic closed to 1 estimates whether  the four HCMC chains converged, values larger than 1.1 suggest that there were convergence issues (Carpenter et al. 2017). 


## Table-S2 legend (tables/Supplementary Tables.xlsx: ciations) 

Parameter coefficients estimate from the effects of the region (temperate vs. tropic) and its interactions with sub-topics in ecology and evolution on the rate of citations. Sub-topics in ecology and evolution are adaptation (AD), climate tolerance (CT), density-dependence (DD), interspecific competition (IC), mimicry (MI), parental care (PC), sexual selection (SS), and speciation (SP). Values of the R statistic closed to 1 estimates whether  the four HCMC chains converged, values larger than 1.1 suggest that there were convergence issues (Carpenter et al. 2017). 


## Table-S3  legend (tables/Supplementary Tables.xlsx: ciations) 

Bayesian model selection testing the importance of the region by sub-topic interaction on the number of publications and on the number of citations. Overall, the models with the region by sub-topic interaction have a better fit to the data, in the number of publications model the interaction has 100% weight of the WAIC, while in the citations model it has 88% of the weight of the WAIC.




 
## References


Carpenter, B., A. Gelman, M. D. Hoffman, D. Lee, B. Goodrich, M. Betancourt, M. Brubaker, J. Guo, P. Li, and A. Riddell. 2017. Stan: A probabilistic programming language. Journal of statistical software 76.


McElreath, R. 2016. Statistical Rethinking: A Bayesian Course with Examples in R and Stan. CRC Press.



@book{mcelreath2016statistical,
  title={Statistical Rethinking: A Bayesian Course with Examples in R and Stan},
  author={McElreath, Richard},
  volume={122},
  year={2016},
  publisher={CRC Press}
}

@article{carpenter2017stan,
  title={Stan: A probabilistic programming language},
  author={Carpenter, Bob and Gelman, Andrew and Hoffman, Matthew D and Lee, Daniel and Goodrich, Ben and Betancourt, Michael and Brubaker, Marcus and Guo, Jiqiang and Li, Peter and Riddell, Allen},
  journal={Journal of statistical software},
  volume={76},
  number={1},
  year={2017},
  publisher={Columbia Univ., New York, NY (United States); Harvard Univ., Cambridge, MA~…}
}

@article{watanabe2010asymptotic,
  title={Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory},
  author={Watanabe, Sumio},
  journal={Journal of Machine Learning Research},
  volume={11},
  number={Dec},
  pages={3571--3594},
  year={2010}
}
