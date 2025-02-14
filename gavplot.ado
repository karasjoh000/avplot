*! create added-variable plots after non-xt estimation
*! version 1.0.0  30apr2024 by John Luke Gallup (jlgallup@pdx.edu)

program define gavplot, rclass sort
	version 11

	local cmd = e(cmd)	// estimation command
	if ("`cmd'"=="") { 
		di as error "no estimates found"
		exit 301
	}
	local ts_cmd "arima arfima arch"
	local sys_cmd "var reg3 sureg mvreg"
	local iv_cmd "ivregress ivreg2"
	local av_cmd "regress `iv_cmd' `ts_cmd' `sys_cmd'"
	if !strpos("`av_cmd'","`cmd'") {
		di as error ///
		   "gavplot not available after command {bf:`cmd'}"
		exit 198
	}
	local multi_eq = strpos("`sys_cmd'","`cmd'")!=0
	
	syntax anything(name=x) [, Level(cilevel) noCI ///
		noCOef CIUnder GENerate(namelist min=2 max=2) ///
		ylim(numlist max=2 ascending) 	///
		xlim(numlist max=2 ascending) 	///
		depvar(varname) ///  only if sys estimator
		Addmeans noDisplay DEBUG *] // undocumented DEBUG option  
			// calculates coefficient from residuals and saves 
			// in r(b_check) to verify calculation of residuals

	local y `e(depvar)'
	if "`x'"!="_cons" fvunab x : `x', max(1)

	_get_gropts, graphopts(`options') gettwoway ///
		getallowed(Rlopts CIOpts addplot)
	local stwowayopts `"`s(twowayopts)'"'
	local gropts `"msize(*.35) pstyle(p1) `s(graphopts)'"'
	local rlopts `"lwidth(medthick) pstyle(p2) `s(rlopts)'"'
	local ciopts `"lpattern(shortdash) pstyle(p2) `s(ciopts)'"'
	local addplot `"`s(addplot)'"'
	_check4gropts rlopts, opt(`rlopts')
	_check4gropts ciopts, opt(`ciopts')

	tempvar touse e_y e_x hat
	tempname b se_b
	preserve

	if "`e(wexp)'" != "" {  // normalize weight variable
		tempvar wvar
		qui gen double `wvar' `e(wexp)'
		if e(wtype) != "fweight" {
			summ `wvar', mean
			qui replace `wvar' = `wvar'/r(mean)
		}
	}
	gen byte `touse' = e(sample)
	
	_ms_lf_info
	if `multi_eq' {
		if "`depvar'"=="" {
			di in red ///
			"depvar() required for multi-equation estimators"
			exit 398
		}
		fvunab depvar :`depvar'
		local eqnum : list posof "`depvar'" in y
		local y `depvar'
		if `eqnum'==0 {
			di in red ///
	"{bf:depvar(`depvar')} is not one of the dependent variables"
			exit 398	 		
		}
		local eqname : word `eqnum' of `e(eqnames)'
		local eqname "`eqname':"
		// check if x in same equation as y
		local X2 = r(varlist`eqnum')
		if r(cons`eqnum') | "`e(method)'"=="2sls" { // has _cons
			tempvar _cons
			gen `_cons' = 1
			if r(cons`eqnum') local X2 "`X2' `_cons'"
		}
		if "`x'"=="_cons" local xtemp `_cons'
		else local xtemp `x'
		local k_x : list posof "`xtemp'" in X2
		if `k_x'==0 {
			di in error /// 
	"variable {bf:`x'} not included in equation {bf:`eqname'}"
			exit 398
		}
		if ("`e(corr)'"=="independent" | "`cmd'"=="mvreg") {
		// process reg3,ols & reg3,2sls & mvreg as single equation
			local X2 : list X2 - x  // remove added var from X2
			local yXvars "`y' `x' `X2'"
			if "`e(method)'"=="2sls" ///
				mata: ivreg("2sls", "`yXvars'",  ///
					"`e(exog)'  `_cons'", "`wvar'",  ///
					"`touse'", "", "`e_y' `e_x'")
			else mata: ols("`yXvars'", "`wvar'", "`touse'",	///
					"`e_y' `e_x'")
		}
		else {  // -reg3- or -sureg-
			local acons = 0
			forvalues e = 1/`r(k_lf)' {
				local Xcount "`Xcount' `r(k`e')'"
				local Xvars "`Xvars' `r(varlist`e')'"
				// add _cons to each equation w/constant
				if (`r(cons`e')') {
					if (!`acons') {
						tempvar _cons
						gen `_cons' = 1
						local acons = 1
					}
					local Xvars "`Xvars' `_cons'"
				}
			}
			if ("`e(endog)'"!="") { // IV
				if (!`acons') {
					tempvar _cons
					gen `_cons' = 1
				}
				local Zvars "`e(exog)' `_cons'"	
			}
			mata: sys("`Xcount'", "`Xvars'", "`Zvars'", ///
					"`wvar'", "`touse'", "`e_y' `e_x'", ///
					`eqnum', `k_x')
		}
	}
	else { // single equation estimator
		if "`depvar'"!="" {
			di in red ///
"{bf:depvar()} option not allowed for single-equation estimators"
			exit 398
		}
		if "`x'"=="`y'" { 
			di in red ///
			"cannot use dependent variable as added variable"
			exit 398
		}
		local X2 = r(varlist)
		if `: list posof "`x'" in X2'==0 {
			di in error /// 
	"variable {bf:`x'} not included in {bf:`cmd'} estimation"
			exit 398
		}
		local X2 : list X2 - x  // remove added var from X2
		if (r(cons1)==1) { // has _cons
			tempvar _cons
			gen `_cons' = 1
			local X2 "`X2' `_cons'"
		}
	}

	scalar `b' 	  = _b[`eqname'`x']
	scalar `se_b' = _se[`eqname'`x']
	
	if ("`cmd'"=="regress") {
		local yXvars "`y' `x' `X2'"
		mata: ols("`yXvars'", "`wvar'", "`touse'", "`e_y' `e_x'")
	}
	else if ("`cmd'"=="ivregress") {
		if ("`e(estimator)'"=="liml") {
			di in red ///
			"{bf:liml} not feasible because of imaginary roots."
			exit 508
		}
		local exogr `e(exogr)'
		local insts `e(insts)'
		local insts : list insts - exogr // remove exogr 
		local X2 "`exogr' `_cons' `e(instd)'"
		local X2 : list X2 - x
		local yXvars "`y' `x' `X2'"
		local Zvars "`exogr' `_cons' `insts'"
		mata: ivreg("`e(estimator)'", "`yXvars'", "`Zvars'", ///
			"`wvar'", "`touse'", "e(kappa)", "`e_y' `e_x'")
	}
	else if ("`cmd'"=="ivreg2") {
		local estim "`e(model)'"
		if ("`estim'"=="ols") ///
			mata: ols("`y' `x' `X2'", "`wvar'", "`touse'", ///
				"`e_y' `e_x'")
		else {
			if ("`estim'"=="liml" | "`estim'"=="kclass") {
				if (e(kclass)>1) {
					di in red ///
			"{bf:`estim'} not feasible because kclass > 1."
					exit 508
				}
				else local estim "kclass"
			}
			else if ("`estim'"=="iv") local estim "2sls"
			local X2 "`e(instd0)' `e(inexog0)' `_cons'"
			local X2 : list X2 - x // remove added var from X2
			local yXvars "`y' `x' `X2'"
			local Zvars "`e(exexog0)' `e(inexog0)' `_cons'"
			mata: ivreg("`estim'", "`yXvars'", "`Zvars'", ///
				"`wvar'", "`touse'", "e(kclass)", "`e_y' `e_x'")
		}
	}
	else if strpos("`ts_cmd'","`cmd'")  /// time series command
		mata: ar("`x'", "`_cons'", "`wvar'", "`touse'", ///
			"`e_y' `e_x'", "`b' `se_b'")
	
	capture local xttl : var label `x'  // capture ∵ x may be fv
	if "`xttl'"=="" local xttl `x'
	capture local yttl : var label `y'
	if "`yttl'"=="" local yttl `y'
	else if "`cmd'"!="regress" local postlbl "*"
	local xttl "e( `xttl'`postlbl' | X`postlbl' )"
	local yttl "e( `yttl'`postlbl' | X`postlbl' )"

	if "`display'"=="" {
		gen `hat' = `b'*`e_x' // regression line

		if ("`ci'"=="") { // create confidence intervals
			tempvar ci_l ci_u			
			tempname t_a
			if ("`cmd'"=="regress") scalar `t_a' = ///
				invttail(`=e(df_r)',(1-`level'/100)/2)
			else scalar `t_a' = -invnormal((1-`level'/100)/2)
			gen `ci_l' = (`b' - `t_a'*`se_b')*`e_x'
			gen `ci_u' = (`b' + `t_a'*`se_b')*`e_x'
		}
		if ("`ylim'`xlim'"!="") { // limit displayed observations
			local yn : word count `ylim'
			if (`yn'>0) {
				if (`yn'==2) {
					local yl1 : word 1 of `ylim'
					local yl2 : word 2 of `ylim'
					local ifyx "`e_y'>`yl1' & `e_y'<`yl2'"
					local ifcil "`ci_l'>`yl1' & `ci_l'<`yl2'"
					local ifciu "`ci_u'>`yl1' & `ci_u'<`yl2'"
					local ifyx2 "`e_y2'>`yl1' & `e_y2'<`yl2'"
				}
				else {
					local ifyx "`e_y'>`ylim'"
					local ifcil "`ci_l'>`ylim'"
					local ifciu "`ci_u'>`ylim'"
					local ifyx2 "`e_y2'>`ylim'"
				}
			}
			local xn : word count `xlim'
			if (`xn'>0) {
				if ("`ifyx'"!="") { // i.e. there are ylims
					local ifyx "`ifyx' & "
					local ifcil "`ifcil' & "
					local ifciu "`ifciu' & "
					local ifyx2 "`ifyx2' & "
				}
				if (`xn'==2) {
					local xlim1 : word 1 of `xlim'
					local xlim2 : word 2 of `xlim'
					local ifyx ///
						"`ifyx' `e_x'>`xlim1' & `e_x'<`xlim2'"
					local ifcil ///
						"`ifcil' `e_x'>`xlim1' & `e_x'<`xlim2'"
					local ifciu ///
						"`ifciu' `e_x'>`xlim1' & `e_x'<`xlim2'"
					local ifyx2 ///
						"`ifyx2' `e_x2'>`xlim1' & `e_x2'<`xlim2'"
				}
				else {
					local ifyx "`ifyx' `e_x'>`xlim'"
					local ifcil "`ifcil' `e_x'>`xlim'"
					local ifciu "`ifciu' `e_x'>`xlim'"
					local ifyx2 "`ifyx2' `e_x2'>`xlim'"
				}
			}
			local ifyx "if `ifyx'"
			local ifcil "if `ifcil'"
			local ifciu "if `ifciu'"
			local ifyx2 "if `ifyx2'"
		}
		if ("`ci'"=="") {	 // after ifyci created	
			local gr_ci "(line `ci_l' `e_x' `ifcil', `ciopts') (line `ci_u' `e_x' `ifciu', `ciopts')"
			if ("`ciunder'"!="") local gr_ci_und `gr_ci'
			else local gr_ci_ovr `gr_ci'
		}
				
		if "`addmeans'"!="" { // add mean values on to e_x and e_y
			quietly {
				tempname ybar xbar
				sum `y' `wgt', meanonly
				scalar `ybar' = r(mean)
				sum `x' `wgt', meanonly
				scalar `xbar' = r(mean)
				replace `e_y' = `e_y' + `ybar'
				replace `hat' = `hat' + `ybar'
				replace `e_x' = `e_x' + `xbar'
				if ("`ci'"=="") {
					replace `ci_l' = `ci_l' + `ybar'
					replace `ci_u' = `ci_u' + `ybar'
				}
			}
			local meanlbl " + mean"
		}
		if ("`coef'" == "") {  // display coefficient estimate in note
			if "`e(vcetype)'"=="Robust" local robust "(robust) "
			if ("`cmd'"=="regress") local t_z "z"
			else local t_z "t"
			local tf : display %5.2f `b'/`se_b'
			local bf : display %9.0g `b'
			local sef : display %9.0g `se_b'
			local note `"note("coef = `bf', `robust'se = `sef', `t_z' = `tf'")"'
		}
		if (`"`addplot'"' != "") local addplot `"(`addplot')"'
		local legend legend(nodraw)
		
		sort `e_x', stable
	
		graph twoway								///
			`gr_ci_und'								///
			(scatter `e_y' `e_x' `ifyx', `gropts') ///
			`scat2'	`gr_ci_ovr'						///
			(line `hat' `e_x' `ifyx', `rlopts')		///
			`addplot'								///
			(, ytitle(`"`yttl'`meanlbl'"') 			///
			  xtitle(`"`xttl'`meanlbl'"')			///
			  `note' `legend' `stwowayopts')
	}
	if "`addmeans'"!="" {
		ret scalar xbar = `xbar'
		ret scalar ybar = `ybar'
	}
	ret scalar se = `se_b'
	ret scalar coef = `b'
	
	if "`debug'"!="" {
		tempname prev_est
		_estimates hold `prev_est', copy
		qui _regress `e_y' `e_x', nocons
		return scalar b_check = _b[`e_x']
		_estimates unhold `prev_est'
	}
	
	if "`generate'"!="" { // exey in mata so protected from -restore-
		mata: exey = st_data(.,"`e_x' `e_y'")
	}	
	restore
	if "`generate'"!="" {
		// check names of generate variables
		capture confirm new variable `generate'
	 	if c(rc) {
			di as err ///
			"variable names in {bf:generate(`generate')} already exist"
			exit 110
		}  
// is "`e(sample)'" causing problems?
		mata: st_store(.,st_addvar("double", tokens("`generate'")), 	       "`e(sample)'", exey)
		mata: mata drop exey
	}
end

mata:
mata set matastrict on

void ols(string scalar yXvars, 	// first X var is the added variable
		string scalar wvar, 	// weight variable if any
		string scalar touse, 	// dummy for included obs
		string scalar evars) // Stata varnames of e_y & e_x
{ // create e_y and e_x from -regress- or -ivreg2, ols- results
	real matrix yX, X2, e
	real colvector w
	real scalar rX
	 
    st_view(yX=., ., yXvars, touse)
	rX = rows(yX)

	if (cols(yX)>2) {  // i.e. there is some X2
		X2 = yX[|1,3\rX,cols(yX)|]
		yX = yX[|1,1\rX,2|] // only y and x1 now
		if (wvar=="") 
		e = yX - quadcross(X2', quadcross(
			invsym(quadcross(X2,X2)), quadcross(X2,yX)))
		else {
			st_view(w=., ., wvar, touse)
			e = yX - quadcross(X2', quadcross(
				invsym(quadcross(X2,w,X2)), quadcross(X2,w,yX)))
		}
	}
	else e = yX		
	st_store(., st_addvar("double", tokens(evars), 1), touse, e)
}

void ivreg(string scalar estim,
			string scalar yXvars, 	// first X is added variable
			string scalar Zvars, 	// instrument variables
			string scalar wvar, 	// weight variable if any
			string scalar touse, 	// dummy for included obs
			string scalar kname, 	// k-class name for liml
			string scalar evars) 	// Stata varnames e_y & e_x
{ // create e_y and e_x from -ivregress- or -ivreg2- results
	real matrix yX, Z, W, D, DZ, E, yX_, X2_, e
	real scalar rX
	real colvector w, d  // weight vector & eigenvector
    st_view(yX, ., tokens(yXvars), touse)
	st_view(Z, ., tokens(Zvars), touse)
	rX = rows(yX)
	
	// calculate (y_tilde,X_tilde) in yX_
	if (wvar=="" & estim=="2sls") {
	   yX_ = quadcross(Z',quadcross(invsym(quadcross(Z,Z)), 
							quadcross(Z,yX)))
	}
	else {
		if (wvar=="") DZ = Z
		else {
			st_view(w, ., wvar, touse)
			DZ = quadcross(diag(w),Z)
		}
		if (estim=="2sls") 
			symeigensystem(quadcross(DZ',
				quadcross(invsym(quadcross(DZ,Z)),DZ')), E, d) 
		else if (estim=="kclass") {
			if (wvar=="") D = I(rows(Z))
			else D = diag(w)
			symeigensystem(D-st_numscalar(kname)*
				(D-quadcross(DZ',quadcross(
					invsym(quadcross(DZ,Z)),DZ'))), E, d)
			if (min(d) < -1e-13) {
				printf(
				`"{err}Added-variables not feasible because "'
					+ `"of imaginary roots.\n"')
				exit(508)
			}
		}
		else { // gmm
			W = st_matrix("e(W)")
			if (st_global("e(cmd)")=="ivreg2") {
				if (cols(DZ)>cols(W)) /// drop fv b. & o. cols
					DZ = select(DZ, colsum(DZ):!=0)
			}
			symeigensystem(quadcross(DZ',quadcross(W,DZ')), E, d) 
		}
		// sqrt, nixing neg. eigvals
		yX_ = quadcross(E',sqrt(d:*(d:>0)),quadcross(E,yX))
	}
	
	if (cols(yX_)>2) {  // i.e. there is some X2
		X2_ = yX_[|1,3\rX,cols(yX_)|]
		yX_ = yX_[|1,1\rX,2|] // only y and x1 now
		e = yX_-quadcross(X2_', quadcross(
			invsym(quadcross(X2_,X2_)), quadcross(X2_,yX_)))
	}
	else e = yX_
	st_store(., st_addvar("double", tokens(evars), 1), touse, e)
}

void sys(string scalar Xcount,	// number of X vars in each equ
		string scalar Xvars,	// all X vars sequenced by equ
		string scalar Zvars,	// all exogenous vars
		string scalar wvar, 	// weight variable if any
		string scalar touse, 	// dummy for included obs
		string scalar evars, 	// Stata varnames e_y & e_x
		real scalar   eqnum,	// eq number of added-variable
		real scalar   k_x)		// order of a.v. in eq's varlist
{
	real matrix X, Z, DZ, D_sqrt, yX_, E
	real rowvector cXf, cXl		// first & last cols of X
	real colvector y, d, e, 
				   rXf, rXl   	// first & last rows of Xd
	real scalar eq, eqtot, n
	
	y = vec(st_data(., st_global("e(depvar)"), touse))
	st_view(X, ., tokens(Xvars), touse)
	n = rows(X)
	if (wvar!="") {
		st_view(w, ., wvar, touse)
		D_sqrt = diag(sqrt(w))
	}
	else D_sqrt = I(n)
	
	if (Zvars!="") {  // if IV 
		st_view(Z, ., tokens(Zvars), touse)
		if (wvar!="") DZ = diag(w)*Z
		else DZ = Z
		X = quadcross(Z',quadcross(invsym(quadcross(DZ,Z)),
			quadcross(DZ,X))) // X becomes Xhat
	}

	eqtot = ustrwordcount(st_global("e(depvar)"))
	rXl = range(1,eqtot,1):*n
	rXf = (0 \ rXl[|1 \ rows(rXl)-1|]):+1
	cXl = runningsum(strtoreal(tokens(Xcount)))
	cXf = (0,cXl[|1 \ cols(cXl)-1|]):+1
	yX_ = J(rows(y), cXl[cols(cXl)], 0)
	for (eq=1; eq<=eqtot; eq++) { // fill in yX_ with pieces of X
		yX_[|rXf[eq],cXf[eq] \ rXl[eq],cXl[eq]|] = 
			X[|1,cXf[eq] \ n,cXl[eq]|]
	}
	yX_ = y,yX_  // prepend y to yX_
	// find eigens of Sigma^-1
	symeigensystem(invsym(st_matrix("e(Sigma)")), E, d)
	// premultiply by sqrt(Sigma^-1), weighted by sqrt(w)
	yX_ = quadcross(
		quadcross(E',sqrt(d:*(d:>0)),E')#D_sqrt, yX_)
	// move a.v. column to second yX_ column
	cx = cXf[eqnum]+k_x
	K = 2..cols(yX_)
	yX_ = yX_[.,(1,cx,select(K, K:!=cx))]
	if (cols(yX_)>2) {  // i.e. there is some X2
		X2_ = yX_[|1,3 \ rows(yX_),cols(yX_)|] // non-av X vars	
		yX_ = yX_[|1,1 \ rows(yX_),2|] // only y and x1 now
		e = yX_-quadcross(X2_', quadcross(
			invsym(quadcross(X2_,X2_)), quadcross(X2_,yX_)))
	}
	else e = yX

	st_store(., st_addvar("double", tokens(evars), 1), touse, 
		e_adj(e, rXf[eqnum], rXl[eqnum]))
}

real matrix e_adj(real matrix e, // unadjusted residuals
					real scalar re1_f, // first eq row
					real scalar re1_l) // last eq row
{  // adjust e_1 values for influence of other sys equations
	real matrix e_1
	real colvector e_x2, e_y2, re, de2
	real rowvector e_1bar
	real scalar n1
	real scalar a_x, a_y
	
	e_1 = e[(re1_f::re1_l),.]
	n1 = rows(e_1)
	e_1bar = mean(e_1)
	re = 1::rows(e)							// row numbers of e
	de2 = re :< re1_f :| re :> re1_l // dummy for e2 obs
	e_x2 = select(e[.,2], de2)
	e_y2 = select(e[.,1], de2)
	a_x = -e_1bar[2] + sqrt(e_1bar[2]^2 + quadcross(e_x2,e_x2)/n1)
	a_y = (quadcross(e_x2,e_y2)/n1 
			- a_x*e_1bar[1])/(a_x + e_1bar[2])
	
	return(e_1 :- (a_y,a_x))
}
				

void ar(string scalar xvar, 		// added variable name
			string scalar cons, 	// constant, if any
			string scalar wvar, 	// weight variable if any
			string scalar touse, 	// dummy for included obs
			string scalar evars, 	// Stata varnames e_y & e_x
			string scalar bse_names) // tempnames for b & s.e.
{ // create variables e_y and e_x from -arima- or -arch- results
	real matrix yXmiss, yX, emiss, e, yX_, S, E, X2_
	real scalar rX, cy
	real colvector nonmiss, wmiss, w, d  // weight vector & eigenvector
	string rowvector ar
"got here 0"
	cmd = st_global("e(cmd)")
"got here 1"  // "trace" messages for debugging Mata
// st_view collect data from stata.
	yX_return = make_yX(xvar, cons, bse_names)
	yX_return
	st_view(yXmiss, ., 
		yX_return, touse)
        // which ones have missing values from lag structure.
	nonmiss = (rowmissing(yXmiss) :== 0)
"got here 1.1"
	st_select(yX, yXmiss, nonmiss)
	if (wvar=="") w = 1
	else {
		st_view(wmiss, ., wvar, touse) // w incl. missing obs
		st_select(w, wmiss, nonmiss)
"got here 1.2"
	}
	rX = rows(yX)
"got here 2"   
    if (cmd == "arima") {
        arima_b_matrix = st_matrix("e(b)")
        arima_ar_max = st_numscalar("e(ar_max)")
        arima_ma_max = st_numscalar("e(ma_max)")
        arima_sigma = st_numscalar("e(sigma)")
        arima_sigma2 = arima_sigma^2 // NOT SURE IF THIS SHOULD BE SQUARED
        // extract ma columns
        arima_ma_matrix = arima_b_matrix[1, (cols(arima_b_matrix) - arima_ma_max)..(cols(arima_b_matrix) - 1)] 
        // append 1 to the start of the matrix
        arima_ma_matrix = J(1,1,1),arima_ma_matrix
        // compute gamma vector
        gamma = J(arima_ma_max + 1, 1, 0)
        for (i = 1; i <= arima_ma_max + 1; i++) {
            arima_ma_right_submatrix = arima_ma_matrix[1, 1..(cols(arima_ma_matrix) - (i - 1))]
            arima_ma_left_submatrix = arima_ma_matrix[1, (1 + (i - 1))..cols(arima_ma_matrix)]
            gamma[i] = arima_sigma2 * (arima_ma_right_submatrix * arima_ma_left_submatrix')
        }
        // append zeros to the end of the gamma vector.
        n_padding = rX - (arima_ma_max + 1)
        v_padding = J(n_padding, 1, 0)
        v_padding
        gamma_padded = gamma \ v_padding
        // compute the covariance matrix using tolpitz
        S = Toeplitz(gamma_padded, gamma_padded')
    } else if (cmd == "arfima") { 
"in arfima code path"
"entering arfima autocov"
        arfima_b_matrix = st_matrix("e(b)")
        arfima_ar_max = st_numscalar("e(ar_max)") 
        arfima_ma_max = st_numscalar("e(ma_max)")
        arfima_ma_matrix = arfima_b_matrix[1, (cols(arfima_b_matrix) - arfima_ma_max-1)..(cols(arfima_b_matrix) - 2)]
        arfima_ma_matrix = J(1,1,1),arfima_ma_matrix
        arfima_sigma2 = arfima_b_matrix[1, cols(arfima_b_matrix)]
        arfima_d = arfima_b_matrix[1, cols(arfima_b_matrix) - 1]
        gammas = arfima_autocov(rX, arfima_sigma2, arfima_d, arfima_ma_max, arfima_ma_matrix)
        S = Toeplitz(gammas', gammas)
        S
        //S = I(rX)
        //S
    } else {
        S = I(rX)
    }
	// calculate (y_tilde,X_tilde) in yX_
	// take sqrt, nixing neg. eigenvals
	symeigensystem(S, E, d) // find eigens of S
	yX_ = quadcross(quadcross(E',sqrt(d:*(d:>0)),E'),w,yX)	
	
	// empty e_y & e_x including missing obs
	st_view(emiss, ., st_addvar("double", tokens(evars)), touse) 
	st_select(e, emiss, nonmiss)	// remove missing rows
	if (cols(yX_)>2) {  // i.e. there is some X2
        // INFO: select submatrix.
		X2_ = yX_[|1,3\rX,cols(yX_)|] // the non-av X variables	
		yX_ = yX_[|1,1\rX,2|] // only y and x1 now
        // X is X1 here.
		e[.,.] = yX_-quadcross(X2_', quadcross(
			invsym(quadcross(X2_,X2_)), quadcross(X2_,yX_)))
	}
	else e[.,.] = yX_
}


real matrix arfima_autocov(real scalar rows, real scalar sigma2, real scalar d, real scalar q, real rowvector theta)
{
    real matrix gamma, phi
    real scalar j, k, sum_phi, sum_gamma, term

    // Initialize gamma and phi vectors
    gamma = J(1, rows, 0)
    phi = J(1, 2*q+1, 0)

    // Compute phi_k
    for (k = -q; k <= q; k++) {
        sum_phi = 0
        for (s = abs(k); s <= q; s++) {
            sum_phi = sum_phi + theta[s+1] * theta[s-abs(k)+1]
        }
        phi[1, k+q+1] = sum_phi
    }

    // Compute gamma_j
    for (j = 0; j < rows; j++) {
        sum_gamma = 0
        for (k = -q; k <= q; k++) {
            term = phi[1, k+q+1] * exp(lngamma(1 - 2*d)) / (exp(lngamma(1 - d))^2)
            term = term * (pochhammer_memo(d, k-j) / pochhammer_memo(1 - d, k-j))
            sum_gamma = sum_gamma + term
        }
        gamma[1, j+1] = sigma2 * sum_gamma
    }

    return(gamma)

}

// Memoized Pochhammer function using an associative cache
real scalar pochhammer_memo(real scalar d, real scalar i)
{
    real matrix cache
	cache
    if (missing(cache)) cache = J(100, 2, .)  // Initialize cache (adjust size as needed)
   
    // Check if value exists in cache
    for (row = 1; row <= rows(cache); row++) {
        if (cache[row, 1] == d & cache[row, 2] == i) {
            return(cache[row, 3])
        }
    }

    // Compute value if not cached
    real scalar result
    if (i == 0) result = 1
    else if (i > 0) result = (d + i - 1) * pochhammer_memo(d, i - 1)
    else result = pochhammer_memo(d, i + 1) / (d + i)

    // Store result in cache
    for (row = 1; row <= rows(cache); row++) {
        if (missing(cache[row, 1])) { // Find empty slot
            cache[row, 1] = d
            cache[row, 2] = i
            cache[row, 3] = result
            break
        }
    }

    return(result)
}



string scalar make_yX(string scalar xvar,	// added variable
						string scalar cons, // _cons varname
						string scalar bse_names) // b & its s.e.
{  // construct varlist of y and X variables for ar()
   //   and extract b and se_b for xvar
	string rowvector Xvars, yXvars, yroot, ars, ebst, bst, bxst
	string scalar yvar, arblock
	real matrix ebse, bse, ar, bx
	real rowvector bse_av
	real scalar arnum

	yvar = st_global("e(depvar)")
	Xvars = tokens(st_global("e(covariates)"))
	if (Xvars=="_NONE") Xvars = J(1,0,"") // void rowvector
	yXvars = J(1,0,"")
	ars = st_global("e(ar)")'
	if (ars!="") {  // add lag operators onto lagged terms
//		ars = subinword(tokens(ars), "1", "") // remove 1 from L1
		arcount = cols(ars)
		for (i=1; i<=arcount; i++) 
			yXvars = yXvars, "L"+ars[i]+"." :+ (yvar,Xvars)
 	}
	else arcount = 0

	yXvars = Xvars, yXvars // add unlagged Xvars
	if (!anyof(yXvars, xvar)) {  // check if av x is in X
		printf(`"{err}Added-variable {bf:%s} not found"', xvar)
		printf("{err} in RHS covariates.\n")
		printf("{err}Make sure L#. is specified if needed.\n")
		exit(322)		
	}
	if (invtokens(yXvars)!=xvar) /// place xvar as the 2nd var
		yXvars = yvar, xvar, select(yXvars, yXvars:!=xvar)
	else yXvars = yvar, xvar
	if (cons!="") yXvars = yXvars, cons
	
	// find estimates b and b_se for xvar
	ebse = st_matrix("e(b)")' // will have columns b and se
	ebse = ebse, sqrt(diagonal(st_matrix("e(V)"))) // add se
	ebst = st_matrixcolstripe("e(b)")
	yroot = ustrsplit(yvar, "[.]") // in case of D.y
	yroot = yroot[cols(yroot)] // keep last part
	if (st_global("e(cmd)")=="arfima") arblock = "ARFIMA"
	else arblock = "ARMA"
	if (arcount) {  // ar terms for y
		ar = select(ebse, ebst[,1]:==arblock)[|1,1 \ arcount,2|] 
		bse = ar
		bst = "L":+ars:+".":+yvar
	}
	else {
		bse = J(0,2,.)
		bst = J(0,1,"")
	}

	// add ARMAX X and lagged X coefficients
	bxst = select(ebst, ebst[,1]:==yroot)[,2]
	if (bxst!="") {  // has _cons or Xvars
		bx = select(ebse, ebst[,1]:==yroot) // X coefs including _cons
		bse = bse \ bx		// Xvars coefficients
		bst = bst \ bxst
		if (Xvars!="" & arcount) { // Xvars ar terms
			if (bxst[rows(bxst)]=="_cons") { // remove _cons row
				bx = bx[|1,1 \ rows(bx)-1,2|] // ∵ not lagged
				bxst = bxst[|1 \ rows(bxst)-1|]
			}
			bxar = J(0,2,.)
			bxarst = J(0,1,"")
			for (i=1; i<=rows(bx); i++) {
				bxar = bxar \ (bx[i,]:*ar)
				bxarst = bxarst \ ("L":+ars:+".":+bxst[i])
			}
			bse = bse \ bxar
			bst = bst \ bxarst
		}
	}

	// match xvar to bst to put estimates b and se_b in Stata
	bse_av = select(bse, xvar:==bst)
	bse_names = tokens(bse_names)
	st_numscalar(bse_names[1], bse_av[1])
	st_numscalar(bse_names[2], bse_av[2])
	return(invtokens(yXvars))	
}

end


