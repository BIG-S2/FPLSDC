# FPLS-DC: Functional partial least squares through distance covariance for imaging genetics. 
This package consists of five functions: fpls_dc and ifpls_dc are the FPLS-DC algorithm and iterative FPLS-DC algorithm, respectively; fplsdc_test and ifplsdc_test are the test procedures corresponding to fpls_dc and ifpls_dc, respectively. To ensure the efficient processing of real data, we have introduced an accelerated version of the fplsdc_test function specifically designed for GWAS (Genome-Wide Association Study). The details of these functions are summarized as follows:

1. fpls_dc: Functional partial least squares through distance covariance(FPLS-DC). FPLS-DC can be used to presents a new feature screening procedure for predictors X.

2. ifpls_dc: Iterative functional partial least squares through distance covariance(IFPLS-DC). IFPLS-DC can be used to presents a variable selection methodology by adding L0 constraint on X. Besides, this function can also update rank of all predictors in the iterative FPLS-DC algorithm and return the index of all ranked predictors X.

3. fplsdc_test: Test of functional partial least squares through distance covariance. The null distribution of statistic is approximated by gamma distribution and the scalar and shape parameter is approximated by permutation procedure. This function returns the p-value vector corresponding to specific predictors when the parameter index.num is determined; otherwise, this function returns the p-value vector corresponding to all predictors.

4. ifplsdc_test: Test of independence by iterative functional partial least squares through distance covariance(IFPLS-DC). The null distribution of statistic is approximated by gamma distribution and the scalar and shape parameter is approximated by permutation procedure. Based on whether the tested variable is in the active set, the procedure of ifplsdc_test is divided into two parts. The specific steps are described in the article. This function returns the p-value vector corresponding to specific predictors when the parameter index.num is determined; otherwise, this function returns the p-value vector corresponding to all predictors.

5. fplsdc_test_gwas: GWAS test of independence for functional partial least squares through distance covariance(FPLS-DC). The null distribution of the statistic is approximated by the gamma distribution, and the scalar and shape parameters are approximated by modifying the referred scalar and shape parameters. These parameters are obtained through a permutation procedure on a referred data set (X1, Y). The probability of X1 being 0, 1, or 2 is 1/3. Y represents the high-dimensional neuroimaging data in the real dataset.

The other details can be found in Wenliang Pan, Yue Shan, Chuang Li, Shuai Huang, Tengfei Li, Yun Li, and Hongtu Zhu, 2023. FPLS-DC: Functional partial least squares through distance covariance for imaging genetics.

