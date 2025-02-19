* test -gavplot-
	
capture set trace off
capture qui log using test_gavplot, replace
//cd "~/Documents/econ/research/stata/ado/personal/gavplot/"
cscript "gavplot for ivregress, etc." adofile gavplot
about
which gavplot

// check whether sort, stable interferes w/ no preserve
// get rid of asterisks in e(x*|X*)


//webuse friedman2, clear
import delimited "data/ma3ar1e0.05.csv", clear
tset t
arima y_t x1 x2, ar(1) ma(1 2 3)

local x_av x1

capture noisily gavplot `x_av', debug nodisplay
quietly {
	if c(rc)==0 {
		noi di as txt "   _b from arima     =  " ///
			as result %18.16f _b[`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`x_av']-r(b_check)) 
		noi assert reldif(_b[`x_av'], r(b_check)) < 1e-1	
	}	
}

local x_av x2

capture noisily gavplot `x_av', debug nodisplay
quietly {
	if c(rc)==0 {
		noi di as txt "   _b from arima     =  " ///
			as result %18.16f _b[`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`x_av']-r(b_check)) 
		noi assert reldif(_b[`x_av'], r(b_check)) < 1e-1	
	}	
}

exit


webuse hsng2, clear
local x_av hsngval

// example from -help reg3-
 webuse klein, clear
// save klein, replace
//  test options ols, 2sls, sure, ireg3, w/ & w/o weights
//use klein, clear
gen taxw = round(taxnetx)
reg3 (wagepriv consump govt capital1) (consump wagepriv wagegovt), 2sls
ereturn list
// coefficient estimates:
matrix list e(b)
set tracedepth 1
// set trace on
local x_av govt
local dv wagepriv
capture noisily gavplot `x_av', depvar(`dv') debug nodisplay
quietly {
	if c(rc)==0 {
		noi di as txt "   _b from reg3     =  " ///
			as result %18.16f _b[`dv':`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`dv':`x_av']-r(b_check)) 
		noi assert reldif(_b[`dv':`x_av'], r(b_check)) < 1e-11	
	}
}
set trace off

exit


/*
// ivregress & ivreg2

webuse hsng2, clear
save hsng2, replace
use hsng2, clear
local x_av hsngval
// model options: iv(or nothing) gmm2s liml fuller(#) kclass(#) coviv(used w/liml or kclass) cue ols(no endog vars)
// ivreg2 rent pcturban hsngval  // ols regression
// ivreg2 rent pcturban (hsngval = faminc i.region)  // 2sls regression
// ivreg2 rent pcturban (hsngval = faminc i.region), liml  // LIML 
// ivreg2 rent pcturban (hsngval = faminc i.region), fuller(11.5) // Fuller
// ivreg2 rent pcturban (hsngval = faminc i.region), kclass(.8)  // kclass
//ivreg2 rent pcturban (hsngval = faminc i.region), cue  // CUE
ivregress gmm rent pcturban (hsngval = faminc i.region) // GMM
ivreg2 rent pcturban (hsngval = faminc i.region), gmm2s  // GMM
ivreg2 rent pcturban (hsngval = faminc i.region), cue
ivreg2 rent pcturban (hsngval = faminc i.region) [pw=pop] // 2sls w/ weights
capture noisily gavplot `x_av', debug nodisplay
quietly {
	if c(rc)==0 {
		noi di as txt "   _b from ivreg2     =  " ///
			as result %18.16f _b[`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`x_av']-r(b_check)) 
		noi assert reldif(_b[`x_av'], r(b_check)) < 1e-13	
	}
}
exit


// -------- arima --------
use friedman2, clear
local x_av m1
arima consump m1 if tin(, 1981q4), arima(1,0,0)
// check _cons & lagged vars w/ ar()

/*
// tests to make:
//  generate + with arima
//  weights
//  check whether e(sample) adapts to 0 weight (make pop[1]=0)
//	factor variables
//  ts variables
//  nocons or hascons
//  check first stage of 2SLS has constant when "`e(constant)'"!=""
//	different variance options
//	collinear variables - dropped variable
//  gmm, igmm
//  any liml w/kappa<=1?
//  ivreg2 options cue, etc.
//  move -ivreg2, ols- to ols()
//  check arima w/ no ar, no ma, xvars, differencing
//  regression constraints?
//  check av for _cons

// example from -help regress-
sysuse auto, clear
regress mpg weight foreign
local x_av weight
capture noisily gavplot `x_av', debug // nodisplay
quietly {
	if c(rc)==0 {
		noi di as txt "   _b from regress     =  " ///
			as result %18.16f _b[`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`x_av']-r(b_check)) 
		noi assert reldif(_b[`x_av'], r(b_check)) < 1e-13	
	}	
}
  
// example from -help arima-
// webuse friedman2, clear
// save friedman2, replace
use friedman2, clear
arima consump m1 if tin(, 1981q4), arima(3,0,0)

// ivregress & ivreg2
webuse hsng2, clear
save hsng2, replace
use hsng2, clear
local x_av "hsngval"
// model options: iv(or nothing) gmm2s liml fuller(#) kclass(#) coviv(used w/liml or kclass) cue OLS(no endog vars)
ivreg2 rent pcturban hsngval `w'  // ols regression
//ivregress 2sls rent pcturban (hsngval = faminc i.region)
ereturn list
set tracedepth 1
//set trace on
	capture noisily gavplot `x_av', debug nodisplay
	quietly {
	if c(rc)==0 {
		noi di as txt "   _b from ivregress     =  " ///
			as result %18.16f _b[`x_av']
		noi di as txt "   _b from reg ey ex =  " ///
			as result %18.16f r(b_check)
		noi di as txt "   diff in _b: " ///
			as result %18.16f (_b[`x_av']-r(b_check)) 
		noi assert reldif(_b[`x_av'], r(b_check)) < 1e-13	
	}
	}
set trace off

exit

// ivregress

webuse hsng2
local x_av "hsngval"
foreach est in 2sls liml gmm {
	foreach w in "" "[w=pop]" {
		ivregress `est' rent pcturban (hsngval = faminc reg1 reg4) `w'
		capture noisily gavplot `x_av', debug nodisplay
		if c(rc)==0 {
			di as txt "   _b from ivregress     =  " as result %18.16f _b[`x_av']
			di as txt "   _b from reg ey ex =  " as result %18.16f r(b_check)
			di as txt "   diff in _b:  " as result %18.16f (_b[`x_av']-r(b_check)) 
//			assert reldif(_b[`x_av'], r(b_check)) < 1e-13	
		}
	}
}
exit

foreach m in fe be re mle pa { 
	di as txt _newline "Model: `m'"

	// x_av not included in preceding -xtreg-
	qui xtreg ln_w age tenure not_smsa, `m'
	capture noisily xtavplot `x_av', debug
	if c(rc)==0 {
		di 
		di as txt "   _b for x not included =  " as result %18.16f r(b_check)
		di as txt "   diff in _b:  " as result %18.16f (scalar(b)-r(b_check)) 
		assert reldif(scalar(b), r(b_check)) < 1e-13	
	}
}

// test for factor variables in varlist
xtreg ln_w c.age#c.age ttl_exp age tenure, fe
xtavplot c.age#c.age

// -------- test that -generate()- option works ----------

local x_av "ttl_exp"
foreach m in be fe re { 
	di as txt _newline "Model: `m'"
	xtavplot `x_av', gen(e_x e_y) debug nodisplay
	noi di as txt "   r(b_check)        =  " as result %18.16f r(b_check)
	_regress e_y e_x, nocons dof(`dof')
	di as txt "   b from xtreg, `m'   =  " as result %18.16f scalar(b)
	di as txt "   _b from e_y on e_x =  " as result %18.16f _b[e_x]
	di as txt "   reldif = " as result %18.16f reldif(_b[e_x],scalar(b)) 
	assert reldif(_b[e_x],scalar(b)) < 1e-13	
	di as txt "   se from xtreg, `m'   =  " as result %18.16f scalar(se)
	di as txt "   _se from e_y on e_x =  " as result %18.16f _se[e_x]
	di as txt "   reldif = " as result %18.16f reldif(_se[e_x],scalar(se)) 
	assert reldif(_se[e_x],scalar(se)) < 1e-13	
	drop e_x e_y
}

// ------ test for varname problems ------

	qui xtreg ln_w age tenure not_smsa, fe
	// no varname
	rcof "xtavplot" == 100

	// too many variables
	rcof "xtavplot `x_av' age" == 103

// ------ test for existence of previous estimates ------

	ereturn clear
	rcof "xtavplot `x_av'" == 301

// ----------- test that INCLUDED x variable is OK -----------

// check for "x was dropped from model"
	gen coll_x = ttl_exp+1
	qui xtreg ln_w coll_x ttl_exp age tenure not_smsa, fe
	if (c(stata_version)<12) rcof "xtavplot coll_x" == 399
	else rcof "xtavplot ttl_exp" == 399
	drop coll_x

	
// ----------- test gavplots -----------

foreach m in fe be re { 
	di as txt _newline "Model: `m'"

	// x_av included in preceding -xtreg-
	qui xtreg ln_w ttl_exp age tenure not_smsa, `m'
	xtavplots
}

// ----------- test -addmeans- option -----------
foreach m in fe be re { 
	di as txt _newline "Model: `m'"

	// x_av included in preceding -xtreg-
	if "`m'"=="fe" local wgt "[aw=birth_yr]"
	else local wgt
	qui xtreg ln_w ttl_exp age tenure not_smsa `wgt', `m'
	capture drop samp
	gen samp = e(sample)
	capture noisily xtavplot `x_av', addmeans
	ret li
	scalar ybar = r(ybar)
	scalar xbar = r(xbar)
	sum ln_w if samp `wgt'
	assert reldif(r(mean),scalar(ybar)) < 1e-13	
	sum `x_av' if samp `wgt'
	assert reldif(r(mean),scalar(xbar)) < 1e-13	
}

display "test_gavplot completed successfully!"

qui log close
window manage forward results
*/


