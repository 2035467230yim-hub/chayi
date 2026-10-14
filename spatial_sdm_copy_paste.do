* Partial dynamic spatial Durbin model, original data, 278 cities, 2009-2024.
* Only lagged structural growth has a spatial Durbin term.
* This file contains ASCII characters only, for legacy Stata editors.
version 18.0
clear all
set more off

* Copy the TWO ORIGINAL DTA files into D:/stata/chayi first.
* Rename copies: raw panel -> data.dta; id/province map -> region.dta.
* Change ROOT if you use another ASCII-only folder.
global ROOT "D:/stata/chayi"
global OUT "$ROOT/spatial_sdm_stru_only"
capture mkdir "$OUT"
capture log close _all
log using "$OUT/run.log", text replace

capture which xsmle
if _rc {
    display as error "Install xsmle first: ssc install xsmle"
    exit 199
}
capture which esttab
if _rc {
    display as error "Install esttab first: ssc install estout"
    exit 199
}

* Build the outcome before selecting cities, exactly as in the baseline file.
use "$ROOT/data.dta", clear
drop if missing(id,year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)
gen double growth_pp=pergdp_real_growth_pct
egen double mean_growth=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-mean_growth)
gen double L1_stru=L.stru
gen double pop_proxy=gdp*10000/pergdp if gdp>0 & pergdp>0
gen double ln_pc_stu=ln(stu/pop_proxy*10000) if stu>0 & pop_proxy>0
gen double pc_retail_10k=retail/pop_proxy if retail>0 & pop_proxy>0
gen double ln_mobile=ln(mobile) if mobile>0
gen double ln_pc_loan=ln(loan/pop_proxy*10000) if loan>0 & pop_proxy>0
merge m:1 id using "$ROOT/region.dta", ///
    keepusing(province) keep(master match) nogen
encode province, gen(province_id)
local controls ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan
egen byte miss_base=rowmiss(abs_gap demo L1_stru `controls' province_id)
gen byte valid=inrange(year,2009,2024) & miss_base==0
bysort id: egen byte nvalid=total(valid)
keep if nvalid==16 & inrange(year,2008,2024)
sort id year
isid id year
xtset id year
assert _N==4726
bysort id: assert _N==17
assert !missing(abs_gap,lon,lat)

* Scale each core variable by its 2009-2024 standard deviation.
quietly summarize demo if inrange(year,2009,2024)
scalar sd_demo=r(sd)
quietly summarize L1_stru if inrange(year,2009,2024)
scalar sd_stru=r(sd)
assert scalar(sd_demo)>0 & scalar(sd_stru)>0
gen double demo_sd=demo/scalar(sd_demo)
gen double stru_sd=L1_stru/scalar(sd_stru)
display as result "SD demo = " %12.6f scalar(sd_demo)
display as result "SD L1_stru = " %12.6f scalar(sd_stru)

* Province-by-year fixed effects; xsmle absorbs city fixed effects.
quietly summarize province_id, meanonly
local np=r(max)
local py ""
forvalues p=1/`np' {
    forvalues t=2010/2024 {
        gen byte py_`p'_`t'=(province_id==`p' & year==`t')
        local py "`py' py_`p'_`t'"
    }
}

* Geographic inverse-square weights from original city coordinates.
preserve
keep if year==2008
sort id
isid id
mkmat lon lat, matrix(COORD)
restore
mata:
    C=st_matrix("COORD")
    n=rows(C)
    W=J(n,n,0)
    pi=acos(-1)
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) if (i!=j) {
            a=sin((C[i,2]-C[j,2])*pi/360)^2+cos(C[i,2]*pi/180)*cos(C[j,2]*pi/180)*sin((C[i,1]-C[j,1])*pi/360)^2
            a=max((0,min((1,a))))
            d=2*6371.0088*asin(sqrt(a))
            if (d<=0) _error(3498)
            W[i,j]=1/d^2
        }
    }
    W=diag(1:/rowsum(W))*W
    st_matrix("W_INV2",W)
end

* Dynamic SDM: lagged outcome, spatial outcome, and spatial lag of stru_sd.
set seed 20261014
xsmle abs_gap demo_sd stru_sd `controls' `py', ///
    wmat(W_INV2) model(sdm) fe type(ind) dlag(1) ///
    durbin(stru_sd) effects vceeffects(sim, nsim(500)) ///
    nolog iterate(500)
count if e(sample)
if r(N)!=4448 {
    display as error "Estimation sample is not 4448; inspect year counts."
    tab year if e(sample)
    exit 459
}
tab year if e(sample)
matrix list e(b)

* Full results: Stata estimates, RTF, CSV, and run log.
estimates store SDM_FULL
estimates save "$OUT/full_estimates.ster", replace
esttab SDM_FULL using "$OUT/full_results.rtf", replace rtf ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N, fmt(0) labels("Estimation sample"))
esttab SDM_FULL using "$OUT/full_results.csv", replace csv ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01)

* Compact results: short- and long-run total effects for both core variables.
local core SR_Total:demo_sd SR_Total:stru_sd ///
    LR_Total:demo_sd LR_Total:stru_sd
esttab SDM_FULL using "$OUT/core_effects.rtf", replace rtf ///
    keep(`core') order(`core') ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    coeflabels(SR_Total:demo_sd "demo short total" ///
               SR_Total:stru_sd "L1_stru short total" ///
               LR_Total:demo_sd "demo long total" ///
               LR_Total:stru_sd "L1_stru long total") ///
    stats(N, fmt(0) labels("Estimation sample"))
esttab SDM_FULL using "$OUT/core_effects.csv", replace csv ///
    keep(`core') order(`core') ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    coeflabels(SR_Total:demo_sd "demo short total" ///
               SR_Total:stru_sd "L1_stru short total" ///
               LR_Total:demo_sd "demo long total" ///
               LR_Total:stru_sd "L1_stru long total")
display as result "Output folder: $OUT"
log close
