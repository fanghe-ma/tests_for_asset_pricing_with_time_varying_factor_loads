# Replication of Tests of asset pricing with time-varying factor loads

## Replication Plan
1. Extend dataset with data from 1963-2025.

2. Implement estimation of factor loadings. Ensure new implementation matches Galvao et al's results when applied on the original dataset. 

3. Implement and validate factor risk premia estimators with 2 different approaches 
    - 3.1: PyTorch penalty based minimization 
    - 3.2: Iterative principle components approach, as described by Galvao et al.

4. Implement estimation of asymptotic variance of factor risk premia estimator, verify via MC and CI coverage

5. Implement calculation of test statistics. Verify asymptotic power and size through MC simulation. In particular 
    - Show trend in power as sample size increases
    - Show trend in size as sample size increases 
    - Show trend in power as heterogeneity increases

6. Implement empirical homogeneity tests 
    - Replicate Galvao et al's empirical test statistics
    - Extend study using extended dataset
