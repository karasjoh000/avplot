discard
import delimited "data/ma3ar1e1.0n50.csv", clear
tset t
arima y_t x1 x2, ar(1) ma(1 2 3)
gavplot x2, debug
return list
