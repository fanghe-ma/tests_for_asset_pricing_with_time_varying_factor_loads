* Galvao, Montes-Rojas, Olmo
* code for data generation
* It requires data_omega_v1_R1.dta
* It creates *.csv matrices to be used in R codes


set more off

* Specify directory where all the m.data is
cd "F:\SDF\realized data full sample"

* Parameters
gl N=47
gl T=206 /*total number of quarters 1963-2014*/

cap program drop running
program define running

qui {
	use data_omega_v1_R1, clear
	gen one=1

	keep if t>=$periodmin & t<=$periodmax
	global T=$periodmax-$periodmin+1

	/*
foreach model in 1 3 5 {
	gen c_model`model'_r=.
	gen c_model`model'_rp=.
	forvalues factor=1(1)`model' {
		gen b_`factor'_model`model'_r=.
		gen b_`factor'_model`model'_rp=.
		}
}
*/
gen id=.

foreach model in 1 3 5 {
gl model=`model'

forvalues id= 1(1)$N {
	replace id=`id' in `id'
		forvalues factor=1(1)`model' {
			reg omega`factor'_`id' L(1/$ar).omega`factor'_`id'
			predict hat
			cap drop Xhat`factor'_`id'
			gen Xhat`factor'_`id'=L.hat
			cap drop hat
			cap drop v`factor'_`id'
			predict v`factor'_`id', resid
			cap drop X`factor'_`id'
			summ L1.omega`factor'_`id'
			gen X`factor'_`id'=L2.omega`factor'_`id'-r(mean)
			cap drop X`factor'_`id'_v
			gen X`factor'_`id'_v=L2.omega`factor'_`id'*L.v`factor'_`id'
		}
		cap drop ones_`id'
		gen ones_`id'=1
}
* ends loop for 47 industries

		* Returns
		gl depvar1=""
		gl depvar2=""
		global regressors=""
		global Xs=""
		global vs=""
		forvalues id= 1(1)$N {
**********************************************
* Modification to substract facors_t*beta_it *
**********************************************
			gl factor_beta=""
			cap drop *_dm
			forvalues f= 1(1)5 {
				summ factor`f'
				gen factor`f'_dm=factor`f'-r(mean)
			}
			global regressors="${regressors} ones_`id'"
			forvalues factor=1(1)`model' {
				global regressors="${regressors} Xhat`factor'_`id'"
				global Xs="${Xs} X`factor'_`id'"
				global vs="${vs} v`factor'_`id'"
				global factor_beta="${factor_beta} - factor`factor'_dm*Xhat`factor'_`id'"
			}
			cap drop depvar2_`id'
			gen depvar2_`id'=aggreturns_`id' + $factor_beta
			gl depvar1="${depvar1} aggreturns_`id'"
			gl depvar2="${depvar2} depvar2_`id'"
		}
		save temp, replace
		drop if Xhat1_1==.
		export delimited $regressors using AX_`model'_$period, replace novarnames
		export delimited $Xs using AXb_`model'_$period, replace novarnames
		export delimited $vs using Av_`model'_$period, replace novarnames
		export delimited $depvar1 using AY1_`model'_$period, replace novarnames
		export delimited $depvar2 using AY2_`model'_$period, replace novarnames
		use temp, replace		
}
* ends model loop

}
* quietly ends here
end



* For model 1 = Sharpe model - 1 factor
* For model 3 = Fama and French (1993) model - 3 factor
* For model 5 = Fama and French (2015) model - 5 factor
gl periodmin=1
gl periodmax=206
gl T=$periodmax-2
gl ar=1
gl period=2014
qui noisily display "**********************************************************************"
qui noisily display "Sample: $periodmin <= t <= $periodmax"
running



