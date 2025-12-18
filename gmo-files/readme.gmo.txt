Antonio Galvao, Gabriel Montes-Rojas and Jose Olmo, "Test of Asset
Pricing with Time-Varying Factor Loadings", Journal of Applied
Econometrics, Vol. 34, No. 5, 2019, pp. 762-778.

All files are stored in the zip file gmo-files.zip.

To replicate the results in the empirical application we need to
follow the following steps:
 
1.- Obtain daily data on 47 Fama-French industry portfolios. We had to
discard two industry portfolios due to data issues.

2.- Use these daily data to construct realized measures of monthly
variances and realized betas.

3.- Do time series regressions (one for each industry portfolio) of
the excess returns on the realized betas.

4.- Do the slope homogeneity tests that we propose in the paper and
the other alternative methods that we consider in the paper.

*****

1.- Daily data on the 47 Fama-French industry portfolios used in the
paper are obtained from Kenneth French's website
(http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)
The dataset relevant for this empirical application is attached under
the name datareturns.txt.

2.- The Matlab file aggfactors.m computes the aggregate factor returns
from high-frequency (daily) data in datareturns.txt

-- The Matlab file datamonthly.m computes the start and end date to
compute quarterly observations from daily observations.  The start and
end date depend on the specific evaluation period, as there can be bank
holidays and weekends that differ across months and years.

-- The Matlab file realizedbetas.m computes the aggregate dynamic betas
from high-frequency (daily) data in datareturns.txt. The output of
this file is the dataset aggreturns.txt that contains the aggregate
quarterly returns and a battery of omegareg*.txt files that contains
the time series of dynamic betas.

3.- GalvaoMontesOlmo_datagen.do uses aggreturns.txt, omegareg*.txt, 
aggfactors.txt to be combined into a single DTA file
data_omega_v1_R1.dta

-- GalvaoMontesOlmo_datamatrix.do runs AR(1) models in eq. (3.7) and
constructs eq. (3.8).  The results are saved into the *.csv files
AY1_*_2014.csv, AX_*_2014.csv, AXb_*_2014.csv, Av_*_2014.csv

4.- GalvaoMontesOlmo.tests.R runs the tests proposed in the paper
using the above *.csv files.

