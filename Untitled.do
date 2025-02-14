discard
import delimited "arima_data.csv", clear
tset t
arfima y_t x1 x2, ar(1) ma(1 2 3)
gavplot x1, debug
return list