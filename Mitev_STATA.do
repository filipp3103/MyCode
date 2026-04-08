
// #1 
// import each file to STATA and conevrt them from .csv to .dta //

import delimited "comtrade_UK_Ger_Fran_Spa_Neth_2020.csv", clear
save "comtrade_UK_Ger_Fran_Spa_Neth_2020.dta", replace 

import delimited "comtrade_UK_Ger_Fran_Spa_Neth_2021.csv", clear
save "comtrade_UK_Ger_Fran_Spa_Neth_2021.dta", replace 

import delimited "comtrade_UK_US_World_2020.csv", clear
save "comtrade_UK_US_World_2020.dta", replace 

import delimited "comtrade_UK_US_World_2021.csv", clear
save "comtrade_UK_US_World_2021.dta", replace 

// now append these files to get the main "dataset" // 
use "comtrade_UK_Ger_Fran_Spa_Neth_2020.dta", clear

append using "comtrade_UK_Ger_Fran_Spa_Neth_2021.dta"
save "dataset1.dta", replace

append using "comtrade_UK_US_World_2020.dta"
save "dataset2.dta", replace

append using "comtrade_UK_US_World_2021.dta"
save "data.dta", replace 

// #2 
// I use a loop to calculate trade valuse of each commodities for each country, so excluding partner "world" //
use "data.dta", clear

levelsof commoditycode, local(commodities)

// I also add a new temporary dataset to store the results //
clear
tempfile results 
save `results', emptyok

foreach i of local commodities {
     use "data.dta"
     quietly summarize tradevalueus if commoditycode == `i' & tradeflowcode == 1 & partnercode != 0
     local sum = r(sum)
     clear
      set obs 1
      gen commoditycode = `i'
      gen importsum = `sum'
   
     append using `results'
      save `results', replace
}

save "summed up.dta", replace

// Change the display format of the variable 'importsum' //
format importsum %15.0f

// Sort and find 5 commodities with the highest value of imports //
gsort -importsum
list commoditycode importsum in 1/5


// #3 
use "data.dta", clear
// to have better looking tables I rename "United States of America" as "USA"
replace partner = "USA" if partner == "United States of America" 
// Imports //
tabstat tradevalueus if tradeflowcode ==1, by(partner) statistics(mean sd min max n)

// Exports //
tabstat tradevalueus if tradeflowcode ==2, by(partner) statistics(mean sd min max n)


// #4
// First, in order to build these graps I find the total value of exports/imports for each country at each period. I summarize trade in all the commodities using loops //

// Imports //

use "data.dta", clear
levelsof period, local(periods)
levelsof partnercode, local(partners)

clear
tempfile results
save `results', emptyok
   
foreach i of local partners {
     foreach j of local periods {
	 use "data.dta"
      quietly summarize tradevalueus if tradeflowcode == 1 & period == `j' & partnercode == `i'

     
	 local sum = r(sum)
      clear
      set obs 1
      gen partner = `i'
      gen period = `j'
      gen importsum = `sum'
      append using `results'
      save `results', replace
}
}

save "results_import.dta", replace

// Exports //

use "data.dta", clear
levelsof period, local(periods)
levelsof partnercode, local(partners)

clear
tempfile results
save `results', emptyok
   
foreach i of local partners {
     foreach j of local periods {
	 use "data.dta"
      quietly summarize tradevalueus if tradeflowcode == 2 & period == `j' & partnercode == `i'

     
	 local sum = r(sum)
      clear
      set obs 1
      gen partner = `i'
      gen period = `j'
      gen exportsum = `sum'
      append using `results'
      save `results', replace
}
}

save "results_export.dta", replace

// Now proceed with graphs // 
// I do a little bit of preparatory work to make graphs look readable

// Convert partner code to the country's name
use "results_import.dta"
tostring partner, generate(partner_)
drop partner
rename partner_ partner

replace partner = "France" if partner == "251"
replace partner = "Spain" if partner == "724"
replace partner = "Germany" if partner == "276"
replace partner = "Netherlands" if partner == "528"
replace partner = "World" if partner == "0"
replace partner = "United States of America" if partner == "842"

// then convert variable "period" to a datetime variable
gen str6 period_str = string(period, "%06.0f")

gen int year = real(substr(period_str, 1, 4))
gen byte month = real(substr(period_str, 5, 2))
gen date_ym = ym(year, month)
format date_ym %tm

// Now do the graph
twoway (line importsum date_ym if partner == "France", lcolor(blue)) (line importsum date_ym if partner == "Spain", lcolor(orange)) (line importsum date_ym if partner == "Netherlands", lcolor(green)) (line importsum date_ym if partner == "Germany", lcolor(red)) (line importsum date_ym if partner == "United States of America", lcolor(black)), legend(order(4 "Germany" 5 "United States of America" 3 "Netherlands" 1 "France" 2 "Spain")) title("UK Imports") ytitle("Imports") xtitle("Date") tlabel(2020m1(3)2021m7, angle(30))
 
 // Now do the same for Exports 
use "results_export.dta", clear
tostring partner, generate(partner_)
drop partner
rename partner_ partner

replace partner = "France" if partner == "251"
replace partner = "Spain" if partner == "724"
replace partner = "Germany" if partner == "276"
replace partner = "Netherlands" if partner == "528"
replace partner = "World" if partner == "0"
replace partner = "United States of America" if partner == "842"

gen str6 period_str = string(period, "%06.0f")
gen int year = real(substr(period_str, 1, 4))
gen byte month = real(substr(period_str, 5, 2))
gen date_ym = ym(year, month)
format date_ym %tm

twoway (line exportsum date_ym if partner == "France", lcolor(blue)) (line exportsum date_ym if partner == "Spain", lcolor(orange)) (line exportsum date_ym if partner == "Netherlands", lcolor(green)) (line exportsum date_ym if partner == "Germany", lcolor(red)) (line exportsum date_ym if partner == "United States of America", lcolor(black)), legend(order(4 "Germany" 5 "United States of America" 3 "Netherlands" 1 "France" 2 "Spain")) title("UK Exports") ytitle("Exports") xtitle("Date") tlabel(2020m1(3)2021m7, angle(30))


// #5 //

// At first let's install necessary additional packages and change dataset "data" a little bit to make it more convenient 

ssc install reghdfe
ssc install ftools

// First, drop extra variables which will not be used in the regressions 
use "data.dta", clear 
drop classification year perioddesc aggregatelevel isleafcode reportercode reporter reporteriso partneriso ndpartnercode ndpartner ndpartneriso modeoftransport modeoftransportcode customs customsproccode qtyunitcode qtyunit qty altqtyunitcode altqtyunit altqty netweightkg grossweightkg ciftradevalueus ciftradevalueus fobtradevalueus flag

// Second, delete all "World" observations
keep if inlist(partner, "France", "Germany", "Netherlands", "Spain", "United States of America")

// Third, log(tradevalues)
gen log_trade = log(tradevalueus)

// Then, add dummy-variables
 gen EU_dummy = 0
 replace EU_dummy = 1 if partner == "France" | partner == "Spain" | partner == "Germany" | partner == "Netherlands"
 
 gen Post_dummy = 0
 replace Post_dummy = 1 if period >= 202101

 gen EUxPost = EU_dummy*Post_dummy
 
 save "data_with_dummies.dta"
 
 
 // And finally, divide the dataset into 2 subsets: Imports and Exports 
 use "data_with_dummies.dta", clear
 keep if inlist(tradeflowcode, 1)
 save "imports_regression_data.dta"
 
  use "data_with_dummies.dta", clear
 keep if inlist(tradeflowcode, 2)
 save "exports_regression_data.dta"
 
 // Run the regressions 
 // Imports 
 use "imports_regression_data.dta"
 reghdfe log_trade EUxPost, absorb(commoditycode partnercode period) cluster(partnercode)
 
 // Exports 
 use "exports_regression_data.dta"
 reghdfe log_trade EUxPost, absorb(commoditycode partnercode period) cluster(partnercode)
