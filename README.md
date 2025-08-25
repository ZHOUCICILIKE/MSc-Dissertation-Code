# MSc-Dissertation-Code
Code for my MSc Data Science dissertation project
# File Descriptions
Here is a description of the key files in this repository:

Gait_MixtureOfMarkovChains.R: The main R script to run the complete analysis. This script loads the data, applies the model, and generates the results presented in the dissertation.

MarkovChainMixtureEM.R: An R script that contains the core functions for implementing the Expectation-Maximization (EM) algorithm for the Mixture of Markov Chains model.

my_mmc_model.rds: A saved R object file containing the final trained model. This allows for quick loading of the model without needing to retrain it.

my_mmc_model.stan: The Stan model file, which defines the statistical model for Bayesian analysis.

mcmc.R: An R script related to the MCMC (Markov Chain Monte Carlo) implementation, likely used with the Stan model.

Helper & Test Scripts:

GenerateMarkovChainMixtureData.R: A script used for generating synthetic data to test the model.

TestMarkovChainMixtureEM.R: A script for testing the implementation of the core EM algorithm.
