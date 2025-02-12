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

# output: 
. gavplot ogap
  got here 1
               make_yX():  3201  vector required
                    ar():     -  function returned error
                 <istmt>:     -  function returned error


```
