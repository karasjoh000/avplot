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
import delimited "arima_data.csv", clear
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

```

# synthetic data generation
```bash 
python arima_gen.py
# creates arima_data.csv file of ar(1) ma(3) process with x1 x2 covariates. 
```
