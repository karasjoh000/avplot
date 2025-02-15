# issues: 
```bash
# in stata shell
cd <git-dir>


# this works: 
use http://www.stata-press.com/data/r12/friedman2, clear
arima consump m2 if tin(, 1981q4), ar(1) ma(1)
gavplot m2


# this does not work
webuse usmacro, clear
arima inflation fedfunds ogap, ar(1) ma(1) 
gavplot ogap
>>  got here 1
>>               make_yX():  3201  vector required
>>                    ar():     -  function returned error
>>                 <istmt>:     -  function returned error


# this does not work
import delimited "data/arima_data.csv", clear
tset t
arima y_t x1 x2, ar(1) ma(3)
gavplot x1, debug
>>  got here 1
>>               make_yX():  3201  vector required
>>                    ar():     -  function returned error
>>                 <istmt>:     -  function returned error


# this works: single covariate
arima y_t x1, ar(1) ma(3)
gavplot x1, debug
return list
>>scalars:
>>            r(b_check) =  1.411650995014425
>>               r(coef) =  1.375391806046776
>>                 r(se) =  .1986195193774165


# make_yX seems to fail when passed into st_view with AR > 1
# note: I debug print the output of make_yX this is where the 
# "y_t x1 x2 ..." output comes from.
arima y_t x1 x2, ar(1 2 3 4) ma(1 2 3 4 5 6 7 8 9 10)
gavplot x1, debug
>>. gavplot x1, debug
>>  got here 0
>>  got here 1
>>  y_t x1 x2 L1 2 3 4.y_t L1 2 3 4.x1 L1 2 3 4.x2 __000007
>>variable L1 not found
>>               st_view():  3598  Stata returned error
>>                    ar():     -  function returned error
>>                 <istmt>:     -  function returned error
>>

# here is the same run with AR = 1
arima y_t x1 x2, ar(1) ma(1 2 3 4 5 6 7 8 9 10)
gavplot x1, debug
>>y_t x1 x2 L1.y_t L1.x1 L1.x2 __000007
```

# synthetic data generation
```bash 
python arima_gen.py
# creates data/ma3ar1.csv file of ar(1) ma(3) process with x1 x2 covariates. 
```
