################################################################################
## RD Designs (Parts I and II)
## NBER Summer Institute Methods Lectures, July 2021
## Author: Matias D. Cattaneo and Rocio Titiunik 
## Python version by Ricardo Masini
## Last update: 12-FEB-2022
################################################################################
## Website: https://rdpackages.github.io/
################################################################################
## RDROBUST: install.packages('rdrobust')
## RDDENSITY: install.packages('rddensity',dependencies=TRUE)
## RDLOCRAND: install.packages('rdlocrand')
################################################################################

import numpy  as np
import pandas  as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from rdrobust import rdbwselect, rdrobust, rdplot
from plotnine import *


################################################################################
## Head Start Data
################################################################################
data = pd.read_csv('headstart.csv')

Y = data['mort_age59_related_postHS']
X = data['povrate60']
Z = data[['census1960_pop', 
       'census1960_pctsch1417',
       'census1960_pctsch534',
       'census1960_pctsch25plus',
       'census1960_pop1417',
       'census1960_pop534',
       'census1960_pop25plus', 
       'census1960_pcturban', 
       'census1960_pctblack']]

C = 59.1984

R = X - C
T = 1*(X>C)

################################################################################
## RDPLOTS
################################################################################
rdplot(Y, X, C, binselect="esmv")
rdplot(Y, X, C, p=1)
rdplot(Y, X, C, nbins=1000)

tempdata = pd.DataFrame(R)
tempdata.rename({'povrate60' : 'v1'}, axis='columns',inplace = True)

plot2 = ggplot(tempdata, aes('v1')) + theme_bw(base_size = 17)
plot2 +=  geom_histogram(tempdata.dropna(), aes(x = 'v1'),
                             breaks = np.arange(min(R),1,1),
                             fill = "blue", color = "black", alpha = 1)
plot2 +=  geom_histogram(tempdata.dropna(), aes(x = 'v1'),
                             breaks = np.arange(0,max(R),1),
                             fill = "red", color = "black", alpha = 1)
plot2 += labs(x = "Score", y = "Number of Observations") 
plot2 += geom_vline(xintercept = 0, color = "black")
plot2

################################################################################
## Replicating Ludwig and Miller (2007, QJE)
################################################################################
## Note: HC0 is not available in Stata, but was used by Ludwig and Miller
ols_data = pd.DataFrame({'Y':Y, 'T':T, 'R':R})
out = smf.ols('Y ~ T + R + R*T',data = ols_data[(-9<=R) & (R<=9)]).fit()
print(out.summary())

out = rdrobust(Y, R, h=9, kernel="uni", vce="hc0"); print(out)
out = rdrobust(Y, X, C, h=9, kernel="uni", vce="hc0"); print(out)

################################################################################
## Local Polynomial Methods
################################################################################
out = rdrobust(Y, X, C); print(out)
out = rdrobust(Y, X, C, h=9, kernel="uni", vce="hc0"); print(out)
out = rdrobust(Y, X, C, h=9, kernel="tri", vce="hc0"); print(out)
out = rdrobust(Y, X, C, h=9, kernel="tri"); print(out)
out = rdrobust(Y, X, C); print(out)

out = rdbwselect(Y, X, C, kernel="uni"); print(out)
out = rdbwselect(Y, X, C, kernel="uni", all=True); print(out)
out = rdbwselect(Y, X, C, all=True); print(out)

out = rdrobust(Y, X, C, bwselect="msetwo"); print(out)

## PLOTS -- LOCAL
out = rdrobust(Y, X, C, h=9, kernel="uni", vce="hc0"); print(out)

out = rdrobust(Y, R, h=9, kernel="uni", vce="hc0"); print(out)

out =  rdrobust(Y, X, C); print(out)
rdplot(Y, R, subset=(-out.bws['left'][0]<= R) & (R <= out.bws['right'][0]),
       binselect="esmv", kernel="triangular",
       h= [out.bws['left'][0], out.bws['right'][0]], p=1)

## RDROBUST with covariates
out  = rdrobust(Y, X, C)
len1 = (out.ci['CI Upper'] - out.ci['CI Lower'])[2]
out = rdrobust(Y, X, C, covs=Z)
len2 = (out.ci['CI Upper'] - out.ci['CI Lower'])[2]
print("CI length change: ", round((len2/len1-1)*100,2), "%")

################################################################################
## Falsification/Validation Methods
################################################################################

## Pre-intervention covariates and placebo outcomes
rdplot(Y, X, C)
rdplot(data['census1960_pop'], X, C)
rdplot(data['census1960_pctsch1417'], X, C)
 
out = rdrobust(Y, X, C, h=9, kernel="uni"); print(out)
out = rdrobust(Y, X, C); print(out)

out = rdrobust(data['census1960_pop'], X, C, h=9, kernel="uni"); print(out)
out = rdrobust(data['census1960_pop'], X, C); print(out)

out = rdrobust(data['census1960_pctsch1417'], X, C, h=9, kernel="uni"); print(out)
out = rdrobust(data['census1960_pctsch1417'], X, C); print(out)

## Placebo cutoff
rdplot(Y, R, p=2, binselect="esmv")
out = rdrobust(Y[R>0], R[R>0], c=3); print(out)

## Recall RD Effect
out = rdrobust(Y, X, C); print(out)
rdplot(Y, R, subset=(-out.bws['left'][0]<= R) & (R <= out.bws['right'][0]),
       binselect="esmv", kernel="triangular",
       h=[out.bws['left'][0], out.bws['right'][0]], p=1)

## Different bandwidths
out = rdbwselect(Y, R, all=True); print(out)
out = rdrobust(Y, R); print(out)
out = rdrobust(Y, R, h=9); print(out)
out = rdrobust(Y, R, h=4); print(out)

## Donut hole
out = rdrobust(Y, R); print(out)
out = rdrobust(Y[(R<=-.1) | (.1 <= R)], R[(R<=-.1) | (.1 <= R)]); print(out)







