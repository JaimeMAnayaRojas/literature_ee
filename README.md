# Analysis of the literature_bias dataset
Here you will find our Bayesian analysis on the number of publications (tropical/temparate) in the web of science (WOS) for nine major topics in ecology and evolution.

# For the methods sections:

We performed a generalised linear mixed model with Poisson distribution on the number of papers (N.papers) published in WOS, as explanatory variables, we used REGION (Tropical vs Temparate) and its interactions with TOPIC (see bellow), and TAXON as random effects in Stan (Carpenter 2017) using the R package Rethingking (Mcelreath 2016). To estimate the parameters,  We used Hamilton Monte Carlo (HMC) chaings (4) and 5000 interaction. All non-adaptive priors used were only weakly informative. HMC chains were verified to be well-mixed and stationary. For model comparasion, we use the Widely Applicable Information Criterion (WAIC). WAIC is  a generilized Bayesian version of AIC and has simmilar interpredation  (Watanabe 2010).

# For the results



# References

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
