* Galvao, Montes-Rojas, Olmo
* code for data generation
* It requires aggreturnsp.txt, aggreturns.txt, omegareg*.txt, aggfactors.txt
* It creates data_omega_v1_R1.dta

set more off

* Specify directory where all the m.data is
*cd ******

* Parameters
gl N=47
gl T=206 /*total number of quarters 1963-2014*/


quietly {

 noisily display "Merging all data into one database..."

 * To get aggreturnsp
local count14=0
forvalues id= 1(1)$N {
	local count=`count14'+2
	local count14=`count'+14
	infix aggreturnsp_`id' `count'-`count14' using aggreturnsp.txt, clear
	gen t=_n
	sort t
	if `id'>1 {
		merge t using auxi
		drop _merge
	}
	sort t
	save auxi, replace
}	



 * To get aggreturns
local count14=0
forvalues id= 1(1)$N {
	local count=`count14'+2
	local count14=`count'+14
	infix aggreturns_`id' `count'-`count14' using aggreturns.txt, clear
	gen t=_n
	sort t
		merge t using auxi
		drop _merge
	sort t
	save auxi, replace
}	

* To get realized covariances, omegas
forvalues id= 1(1)$N {
	infix omega1_`id' 2-16 omega2_`id' 18-32 omega3_`id' 34-48 omega4_`id' 50-64 omega5_`id' 66-80 using omegareg`id'.txt, clear
	gen t=_n
	sort t
	merge t using auxi
	drop _merge
	sort t
	save auxi, replace
}	


* To get factors, factors
	infix factor1 2-16 factor2 18-32 factor3 34-48 factor4 50-64 factor5 66-80 using aggfactors.txt, clear
	gen t=_n
	sort t
	merge t using auxi
	drop _merge
	sort t
	save auxi, replace


tsset t
save data_omega_v1_R1, replace
}
