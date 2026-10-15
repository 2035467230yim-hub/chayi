* ASCII-only Stata 18 file: one new moderator and one candidate threshold.
* The first two moderator models and retail threshold are in your earlier files.
* This file uses the same 278 cities, 2009-2024, city/year FE and five
* controls as your uploaded moderation and threshold tables.
* Smaller abs_gap means a smaller regional growth-rate disparity.
version 18.0
clear all
set more off

* Original source folder. Keep the Unicode escapes if the folder has not moved.
local root_esc "D:/Desktop/\u5de5\u4f5c/\u7814\u7a76\u751f\u5e08\u95e8\u5de5\u4f5c/\u4e13\u5229\u8f6c\u79fb/\u3010\u6700\u7ec8\u6210\u529f\u3011\u4e13\u5229\u8f6c\u79fb\u4e0e\u533a\u57df\u7ecf\u6d4e\u589e\u957f\u30102009-2024\u5e74\u3011"
local root=ustrunescape("`root_esc'")
local raw_file=ustrunescape("\u5220\u9664\u90e8\u5206\u57ce\u5e02\u548c\u7f29\u5c3e1%\u7684\u6570\u636e.dta")
capture mkdir "D:/stata"
global OUT "D:/stata/chayi_mechanism_add"
capture mkdir "$OUT"
capture log close _all
log using "$OUT/full_results.log", text replace
capture confirm file "`root'/`raw_file'"
if _rc {
    display as error "Original DTA not found. Edit root_esc at the top."
    exit 601
}
foreach cmd in reghdfe xthreg esttab {
    capture which `cmd'
    if _rc {
        display as error "Missing command: `cmd'"
        exit 199
    }
}

* Build the outcome and one-year lags before restricting the panel.
use "`root'/`raw_file'", clear
drop if missing(id,year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)
gen double growth_pp=pergdp_real_growth_pct
egen double annual_mean=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-annual_mean)
gen double L1_stru=L.stru
gen double L1_gdppc_real=L.pergdp_real
gen double pop_proxy=gdp*10000/pergdp if gdp>0 & pergdp>0
gen double ln_pc_stu=ln(stu/pop_proxy*10000) if stu>0 & pop_proxy>0
gen double ln_pc_lib=ln(lib/pop_proxy*10000) if lib>0 & pop_proxy>0
gen double ln_pc_hosp=ln(hosp/pop_proxy*10000) if hosp>0 & pop_proxy>0
gen double ln_retail=ln(retail) if retail>0
gen double ln_pc_mobile=ln(mobile/pop_proxy*10000) if mobile>0 & pop_proxy>0
gen double road_pc=road/pop_proxy if road>=0 & pop_proxy>0
gen double L1_gdppc_10k=L1_gdppc_real/10000

local controls ln_pc_stu ln_pc_lib ln_pc_hosp ln_retail ln_pc_mobile
egen byte miss_main=rowmiss(abs_gap demo L1_stru `controls')
gen byte valid=inrange(year,2009,2024) & miss_main==0
bysort id: egen byte nvalid=total(valid)
keep if valid & nvalid==16
sort id year
isid id year
xtset id year
assert _N==4448
bysort id: assert _N==16
assert !missing(road_pc,L1_gdppc_10k)

* MODERATION 3: road mileage per population, an accessibility proxy.
* Shortage = mean minus standardized road_pc. Greater shortage means fewer
* road kilometers per resident. Interaction p values test slope differences.
quietly summarize road_pc
scalar road_mean=r(mean)
scalar road_sd=r(sd)
assert scalar(road_sd)>0
gen double road_short=(scalar(road_mean)-road_pc)/scalar(road_sd)
gen double demo_road=demo*road_short
gen double stru_road=L1_stru*road_short
reghdfe abs_gap demo L1_stru road_short demo_road stru_road ///
    `controls', absorb(id year) vce(cluster id)
test demo_road
local p_demo_inter=r(p)
test stru_road
local p_stru_inter=r(p)
estadd scalar P_demo_inter=`p_demo_inter'
estadd scalar P_stru_inter=`p_stru_inter'
estimates store ROAD_MOD
quietly _pctile road_short, p(25 75)
local short25=r(r1)
local short75=r(r2)
lincom L1_stru+`short25'*stru_road
lincom L1_stru+`short75'*stru_road
esttab ROAD_MOD using "$OUT/road_moderation_full.rtf", replace rtf ///
    b(6) se(6) star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N P_demo_inter P_stru_inter, fmt(0 4 4))
esttab ROAD_MOD using "$OUT/road_moderation_core.csv", replace csv ///
    keep(demo L1_stru road_short demo_road stru_road) ///
    b(6) se(6) star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N P_demo_inter P_stru_inter, fmt(0 4 4))

* THRESHOLD 2 CANDIDATE: previous-year real per-capita GDP level.
* Only L1_stru changes slope; demo stays in the common linear component.
* xthreg's own bootstrap p value is the decision rule for a threshold.
quietly tabulate year, generate(_yr)
local yearfe ""
forvalues j=2/16 {
    local yearfe "`yearfe' _yr`j'"
}
xthreg abs_gap demo ln_retail ln_pc_stu ln_pc_lib ln_pc_hosp ///
    ln_pc_mobile L1_gdppc_10k `yearfe', ///
    rx(L1_stru) qx(L1_gdppc_10k) ///
    thnum(1) trim(.05) grid(200) bs(1000) ///
    nobslog vce(cluster id)
matrix TH=e(Thrss)
matrix F=e(Fstat)
local cut=TH[1,1]
local ci_lo=TH[1,4]
local ci_hi=TH[1,5]
local p_boot=F[1,4]
display as result "Lagged GDP threshold = " %9.5f `cut'
display as result "xthreg bootstrap threshold p = " %7.4f `p_boot'

* Conditional slopes after the selected cut; their p values are descriptive.
gen byte high=(L1_gdppc_10k>`cut')
gen double stru_low=L1_stru*(high==0)
gen double stru_high=L1_stru*(high==1)
xtreg abs_gap demo ln_retail ln_pc_stu ln_pc_lib ln_pc_hosp ///
    ln_pc_mobile L1_gdppc_10k ///
    stru_low stru_high i.year, fe vce(cluster id)
quietly test stru_low=stru_high
local p_stru_diff=r(p)
estadd scalar Cut=`cut'
estadd scalar Bootstrap_p=`p_boot'
estadd scalar P_stru_diff=`p_stru_diff'
estimates store GDP_THRESHOLD
esttab GDP_THRESHOLD using "$OUT/gdp_threshold_full.rtf", replace rtf ///
    b(6) se(6) star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N Cut Bootstrap_p P_stru_diff, fmt(0 5 4 4))
esttab GDP_THRESHOLD using "$OUT/gdp_threshold_core.csv", replace csv ///
    keep(demo stru_low stru_high) ///
    b(6) se(6) star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N Cut Bootstrap_p P_stru_diff, fmt(0 5 4 4))
display as result "Output folder: $OUT"
log close
