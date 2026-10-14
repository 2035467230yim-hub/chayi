* ASCII-only Stata 18 do-file. Chinese source paths are decoded at runtime.
* No renaming or copying of the original DTA files is required.
* Partial dynamic SDM; 278 cities, 2009-2024; only stru has W*x.
version 18.0
clear all
set more off

* The next three names decode to the same Chinese paths as your baseline do-file.
* Edit root_esc if your data folder has moved. Keep the \u escapes.
local root_esc "D:/Desktop/\u5de5\u4f5c/\u7814\u7a76\u751f\u5e08\u95e8\u5de5\u4f5c/\u4e13\u5229\u8f6c\u79fb/\u3010\u6700\u7ec8\u6210\u529f\u3011\u4e13\u5229\u8f6c\u79fb\u4e0e\u533a\u57df\u7ecf\u6d4e\u589e\u957f\u30102009-2024\u5e74\u3011"
local root = ustrunescape("`root_esc'")
global ROOT "`root'"
local raw_file = ustrunescape("\u5220\u9664\u90e8\u5206\u57ce\u5e02\u548c\u7f29\u5c3e1%\u7684\u6570\u636e.dta")
local map_file = ustrunescape("\u7701\u4efd\u57ce\u5e02id\u5bf9\u7167.dta")
* Keep results in an ASCII-only folder, even if source paths contain Chinese.
capture mkdir "D:/stata"
global OUT "D:/stata/chayi_results"
capture mkdir "$OUT"
capture log close _all
log using "$OUT/run.log", text replace

capture confirm file "$ROOT/`raw_file'"
if _rc {
    display as error "Original panel data not found at the decoded source path."
    exit 601
}
capture confirm file "$ROOT/`map_file'"
if _rc {
    display as error "Original city/province map not found at the decoded source path."
    exit 601
}

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
use "$ROOT/`raw_file'", clear
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
merge m:1 id using "$ROOT/`map_file'", ///
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
* First-run setting: 100 MC draws for effect SEs; use 500 for final reporting.
set seed 20261014
display as result "Starting xsmle; the effect-SE simulation can take several minutes."
xsmle abs_gap demo_sd stru_sd `controls' `py', ///
    wmat(W_INV2) model(sdm) fe type(ind) dlag(1) ///
    durbin(stru_sd) effects vceeffects(sim, nsim(100)) ///
    iterate(500)
display as result "xsmle and simulated effect SEs finished."
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
