*************************************************************************
** RD Designs (Parts I and II)
** NBER Summer Institute Methods Lectures, July 2021
** Author: Matias D. Cattaneo and Rocio Titiunik
** Last update: 26-JUL-2021
*************************************************************************
** Website: https://rdpackages.github.io/
*************************************************************************
** RDROBUST: net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
** RDDENSITY: net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
** RDLOCRAND: net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
** RDMULTI: net install rdmulti, from(https://raw.githubusercontent.com/rdpackages/rdmulti/master/stata) replace
** RDPOWER: net install rdpower, from(https://raw.githubusercontent.com/rdpackages/rdpower/master/stata) replace
*************************************************************************
** LPDENSITY: net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace
*************************************************************************
clear all
set more off
set linesize 200

*************************************************************************
** Empirical Illustration: Ludwig and Miller (2007, QJE)
*************************************************************************
use headstart, clear
des
gl Y mort_age59_related_postHS
gl X povrate60
gl Z census1960_pop census1960_pctsch1417 census1960_pctsch534 ///
     census1960_pctsch25plus census1960_pop1417 census1960_pop534 ///
	 census1960_pop25plus census1960_pcturban census1960_pctblack

gl C 59.1984

gen R = $X - $C

*************************************************************************
** Summary Statistics
*************************************************************************
sum $Y $X mort_age25plus_related_postHS $Z
gen T = ($X >= $C)
ttest $Y, by(T)

*************************************************************************
** RDPLOTS
*************************************************************************
rdplot $Y $X, c($C) binselect(esmv) graph_options(name("Default"))
rdplot $Y $X, c($C) p(1) graph_options(name("RDPLOT_p1"))
rdplot $Y $X, c($C) nbins(1000)

twoway (histogram $X if $X < $C, freq width(3) bcolor(red)) ///
       (histogram $X if $X >= $C, freq width(3) bcolor(blue) xline($C))

	   
*************************************************************************
** Replicating Ludwig and Miller (2007, QJE)
*************************************************************************
** Note: HC0 is not available in Stata, but was used by Ludwig and Miller
reg $Y T##c.R if abs(R)<=9, robust

rdrobust $Y R, h(9) kernel(uni) vce(hc0)
rdrobust $Y $X, c($C) h(9) kernel(uni) vce(hc0)

rdrandinf $Y R, c(0) wl(-9) wr(9) p(1) reps(10)
rdrandinf $Y $X, c($C) wl(50.1984) wr(68.1984) p(1) reps(10)

*************************************************************************
** Local Randomization Methods
*************************************************************************
rdrandinf $Y $X, c($C) wl(57) wr(61) reps(1000)

rdrandinf $Y $X, c($C) wl(57) wr(61) reps(1000) stat(all)

rdwinselect $X census1960_pop census1960_pcturban census1960_pctblack, ///
            c(59.1984) obsmin(10) wobs(5) nwindows(20) seed(43)

rdplot $Y $X if 57<=abs($X) & abs($X)<=61, c($C) nbins(1000) p(0) graph_options(name("RDPLOT_LR"))



*************************************************************************
** Local Polynomial Methods
*************************************************************************
rdrobust $Y $X, c($C)
rdrobust $Y $X, c($C) h(9) kernel(uni) vce(hc0)
rdrobust $Y $X, c($C) h(9) kernel(tri) vce(hc0)
rdrobust $Y $X, c($C) h(9) kernel(tri)
rdrobust $Y $X, c($C)

rdbwselect $Y $X, c($C) kernel(uni)
rdbwselect $Y $X, c($C) kernel(uni) all
rdbwselect $Y $X, c($C) all

rdrobust $Y $X, c($C) bwselect(msetwo)

** PLOTS -- LOCAL
rdrobust $Y $X, c($C) h(9) kernel(uni) vce(hc0)

rdrobust $Y R, h(9) kernel(uni) vce(hc0)
rdplot $Y R if -e(h_l) <= R & R <= e(h_r), ///
       binselect(esmv) kernel(uniform) h(`e(h_l)' `e(h_r)') p(1) nbins(1000) ///
       graph_options(title("RD Plot") ///
                     ytitle(Child mortality) ///
                     xtitle(Poverty index) ///
                     graphregion(color(white)) ///
					 name("RDPlotLocal1000"))

** RDROBUST with covariates
rdrobust $Y $X, c($C)
local len = `e(ci_r_rb)' - `e(ci_l_rb)'
rdrobust $Y $X, c($C) covs($Z)
display "CI length reduction: " round(((`e(ci_r_rb)'-`e(ci_l_rb)')/`len'-1)*100,.01) "%"

*************************************************************************
** Falsification/Validation Methods
*************************************************************************
graph drop _all

** Density Tests: Binomial and Continuity-Based
twoway (histogram $X if $X < $C, freq width(2) bcolor(red)) ///
       (histogram $X if $X >= $C, freq width(2) bcolor(blue) xline($C)), ///
	   name("Histogram") leg(off)

rddensity $X, c($C) plot graph_opt(name("DensityTest") leg(off))

rddensity $X, c($C)

rdwinselect $X, c($C)

tab T if abs($X-$C) <= 1
bitesti 65 35 1/2

** Pre-intervention covariates and placebo outcomes
rdplot $Y $X, c($C) graph_options(name("Outcome"))
rdplot census1960_pop $X, c($C) graph_options(name("PretreatCovariate"))
rdplot mort_age25plus_related_postHS $X, c($C) graph_options(name("PlaceboOutcome"))
 
rdrobust $Y $X, c($C) h(9) kernel(uni)
rdrobust census1960_pop $X, c($C) h(9) kernel(uni)
rdrobust mort_age25plus_related_postHS $X, c($C) h(9) kernel(uni) vce(hc0)

** Placebo cutoff
rdplot $Y R, p(2) binselect(esmv)
rdrobust $Y R if R>0, c(2) kernel(uni) vce(hc0)
rdrobust $Y R if R>0, c(3) kernel(uni) vce(hc0)

** Recall RD Effect
rdrobust $Y R, h(9) kernel(uni) vce(hc0)
rdplot $Y R if -e(h_l) <= R & R <= e(h_r),  ///
       binselect(esmv) kernel(uniform) h(`e(h_l)' `e(h_r)') p(1) ///
       graph_options(title("RD Plot") ///
                     ytitle(Child mortality) ///
                     xtitle(Poverty index) ///
                     graphregion(color(white)) ///
					 name("RDPlotLocal"))

** Different bandwidths
rdbwselect $Y R, kernel(uni) all
rdrobust $Y R, h(9) kernel(uni) vce(hc0)
rdrobust $Y R, h(4) kernel(uni) vce(hc0)

** Donut hole
rdrobust $Y R, h(9) kernel(uni) vce(hc0)
rdrobust $Y R if abs(R)>=0.10, h(9) kernel(uni) vce(hc0)





