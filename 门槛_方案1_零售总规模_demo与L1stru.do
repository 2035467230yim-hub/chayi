* Scheme 1: retail-market scale as q; both demo and L1_stru vary across one threshold.
version 18.0
clear all
set more off
set maxvar 32767
global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global OUTDIR "C:\Users\Nfoh\Documents\ChatGPT\区域经济差异"
capture log close _all
log using "$OUTDIR\方案1_零售总规模.log", text replace
foreach cmd in xthreg esttab {
    capture which `cmd'
    if _rc {
        display as error "缺少Stata命令：`cmd'"
        log close
        exit 199
    }
}
use "$DATA", clear
foreach v in id year pergdp pergdp_real_growth_pct demo stru gdp ///
    stu lib hosp mobile inv retail {
    capture confirm variable `v'
    if _rc {
        display as error "原始数据缺少变量：`v'"
        log close
        exit 111
    }
}
drop if missing(id) | missing(year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)
gen double growth_pp=pergdp_real_growth_pct
gen double L1_stru=L.stru
gen double pop_proxy=gdp*10000/pergdp if gdp>0 & pergdp>0
gen double ln_pc_stu=ln(stu/pop_proxy*10000) if stu>0 & pop_proxy>0
gen double ln_pc_lib=ln(lib/pop_proxy*10000) if lib>0 & pop_proxy>0
gen double ln_pc_hosp=ln(hosp/pop_proxy*10000) if hosp>0 & pop_proxy>0
gen double ln_retail=ln(retail) if retail>0
gen double ln_pc_mobile=ln(mobile/pop_proxy*10000) if mobile>0 & pop_proxy>0
gen double inv_pc=inv/pop_proxy*10000 if inv>=0 & pop_proxy>0
gen double ln_inv_pc=ln(1+inv_pc) if inv_pc>=0
gen double ln_gdp=ln(gdp) if gdp>0
global CTRL ln_pc_stu ln_pc_lib ln_pc_hosp ln_pc_mobile
egen byte _miss=rowmiss(growth_pp demo L1_stru ln_retail $CTRL ln_inv_pc ln_gdp)
gen byte _valid=inrange(year,2009,2024) & _miss==0
bysort id: egen int _ny=total(_valid)
keep if _valid & _ny==16
drop _miss _valid _ny
egen double _mean=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-_mean)
quietly count
local N=r(N)
quietly count if year==2009
local G=r(N)
display as result "样本城市=`G'；观测值=`N'"
quietly tabulate year, generate(_yr)
global YEARFE ""
forvalues j=2/16 {
    global YEARFE "$YEARFE _yr`j'"
}
quietly xthreg abs_gap ln_retail $CTRL $YEARFE, ///
    rx(demo L1_stru) qx(ln_retail) thnum(1) trim(.05) grid(200) bs(1000) ///
    nobslog vce(cluster id)
matrix TH=e(Thrss)
matrix F=e(Fstat)
local th=TH[1,1]
local lo=TH[1,4]
local hi=TH[1,5]
local pb=F[1,4]
gen byte _high=(ln_retail>`th')
gen double _demo_lo=demo*(_high==0)
gen double _demo_hi=demo*(_high==1)
gen double _stru_lo=L1_stru*(_high==0)
gen double _stru_hi=L1_stru*(_high==1)
quietly xtreg abs_gap ln_retail _demo_lo _demo_hi _stru_lo _stru_hi $CTRL i.year, fe vce(cluster id)
quietly test _demo_lo=_demo_hi
local pd=r(p)
quietly test _stru_lo=_stru_hi
local ps=r(p)
estadd scalar Threshold=`th'
estadd scalar CI_Lower=`lo'
estadd scalar CI_Upper=`hi'
estadd scalar Bootstrap_p=`pb'
estadd scalar Demo_diff_p=`pd'
estadd scalar Stru_diff_p=`ps'
display as result "零售总规模门槛=`th'; bootstrap p=`pb'; demo slope-diff p=`pd'; L1_stru slope-diff p=`ps'"
esttab using "$OUTDIR\结果表_方案1_零售总规模.rtf", replace rtf ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(_demo_lo _demo_hi _stru_lo _stru_hi) order(_demo_lo _demo_hi _stru_lo _stru_hi) ///
    coeflabels(_demo_lo "demo（门槛下）" _demo_hi "demo（门槛上）" ///
        _stru_lo "L1_stru（门槛下）" _stru_hi "L1_stru（门槛上）") ///
    mtitles("零售总规模") ///
    stats(N Threshold CI_Lower CI_Upper Bootstrap_p Demo_diff_p Stru_diff_p r2_w, ///
        fmt(0 4 4 4 4 4 4 4) labels("观测值" "门槛值" "95%区间下界" "95%区间上界" ///
        "1000次bootstrap门槛p值" "demo两侧斜率差异p值" "L1_stru两侧斜率差异p值" "组内R²")) ///
    title("门槛效应检验：零售总规模") ///
    addnotes("因变量为实际人均GDP增长率相对年度均值的绝对离差；demo与L1_stru均允许跨门槛改变斜率。", ///
        "控制人力资本、图书馆、医疗、人均移动电话；城市与年份固定效应，标准误按城市聚类。", ///
        "门槛存在依据1000次bootstrap检验；组间斜率差异检验为补充。")
log close
