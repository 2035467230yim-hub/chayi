* ============================================================================
* 安慰剂检验：年份内随机置换、城市内随机年份错配；各1000次
* 汇报随机回归 t 统计量分布，并以真实基准 t 值作竖线。
* 设定严格匹配最终基准：abs_gap，demo、L1_stru，四控制变量，城市FE+省份×年份FE，城市聚类。
* 只输出 do 对应的 log、CSV、PNG；不保存持久 DTA。
* ============================================================================
version 18.0
clear all
set more off
set maxvar 32767

global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global REGIONMAP "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\省份城市id对照.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\反+安 最后"
global REPS 1000

capture mkdir "$OUTDIR"
capture log close _all
log using "$OUTDIR\安慰剂检验_1000次_t值分布.log", text replace

foreach cmd in reghdfe {
    capture which `cmd'
    if _rc {
        display as error "缺少Stata命令：`cmd'。请先安装 reghdfe。"
        log close
        exit 199
    }
}
capture confirm file "$DATA"
if _rc {
    display as error "找不到面板数据：$DATA"
    log close
    exit 601
}
capture confirm file "$REGIONMAP"
if _rc {
    display as error "找不到省份城市对照表：$REGIONMAP"
    log close
    exit 601
}

use "$DATA", clear
drop if missing(id,year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)

* 因变量、核心解释变量和最终基准四项控制变量。
gen double growth_pp=pergdp_real_growth_pct
egen double mean_growth=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-mean_growth)
gen double L1_stru=L.stru
gen double pop_proxy=gdp*10000/pergdp if gdp>0 & pergdp>0
gen double pc_stu=stu/pop_proxy*10000 if stu>0 & pop_proxy>0
gen double ln_pc_stu=ln(pc_stu) if pc_stu>0
gen double pc_retail=retail/pop_proxy*10000 if retail>0 & pop_proxy>0
gen double pc_retail_10k=pc_retail/10000 if pc_retail>=0
gen double ln_mobile=ln(mobile) if mobile>0
gen double pc_loan=loan/pop_proxy*10000 if loan>0 & pop_proxy>0
gen double ln_pc_loan=ln(pc_loan) if pc_loan>0

merge m:1 id using "$REGIONMAP", keepusing(province) keep(master match) nogen
encode province, gen(province_id)
egen long province_year=group(province_id year)
global CTRL "ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan"

* 与最终基准锁定同一共同样本：278城市×2009—2024年=4448条。
egen byte miss_base=rowmiss(abs_gap demo L1_stru $CTRL province_id)
gen byte base_sample=inrange(year,2009,2024) & miss_base==0
bysort id: egen int n_base=total(base_sample)
keep if n_base==16 & base_sample==1
quietly count
if r(N)!=4448 {
    display as error "共同样本应为4448条，当前为" r(N) "；停止以免与基准错配。"
    log close
    exit 459
}
quietly levelsof id, local(cityids)
local ncity : word count `cityids'
if `ncity'!=278 {
    display as error "共同样本应为278个城市，当前为 `ncity'；停止以免与基准错配。"
    log close
    exit 459
}
display as result "安慰剂共同样本：`ncity'个城市，4448条观测。"

* 真实基准系数、标准误和t值。
quietly reghdfe abs_gap demo L1_stru $CTRL, ///
    absorb(id province_year) vce(cluster id)
scalar b_demo_real=_b[demo]
scalar se_demo_real=_se[demo]
scalar t_demo_real=_b[demo]/_se[demo]
scalar p_demo_real=2*ttail(e(df_r),abs(scalar(t_demo_real)))
scalar b_stru_real=_b[L1_stru]
scalar se_stru_real=_se[L1_stru]
scalar t_stru_real=_b[L1_stru]/_se[L1_stru]
scalar p_stru_real=2*ttail(e(df_r),abs(scalar(t_stru_real)))
display as result "真实基准：demo b=" %10.7f scalar(b_demo_real) ///
    " t=" %8.4f scalar(t_demo_real) " p=" %7.4f scalar(p_demo_real)
display as result "真实基准：L1_stru b=" %10.7f scalar(b_stru_real) ///
    " t=" %8.4f scalar(t_stru_real) " p=" %7.4f scalar(p_stru_real)

tempfile base_panel draws_year draws_city summary
save `base_panel', replace
tempname hsum
postfile `hsum' str32 placebo int repetitions ///
    double real_b_demo real_t_demo p_demo_t_two p_demo_t_left ///
    real_b_stru real_t_stru p_stru_t_two p_stru_t_left ///
    both_negative_sig10_share using `summary', replace

* ---------------------- 方案1：在同一年内随机错配城市 ----------------------
tempname hyear
postfile `hyear' int rep ///
    double b_demo se_demo t_demo p_demo b_stru se_stru t_stru p_stru ///
    byte both_neg10 using `draws_year', replace
set seed 20261009
forvalues r=1/$REPS {
    preserve
        use `base_panel', clear
        tempfile main draw
        tempvar u rank u2
        gen double `u'=runiform()
        bysort year (`u'): gen long `rank'=_n
        save `main', replace

        keep year `rank' demo L1_stru
        rename demo placebo_demo
        rename L1_stru placebo_L1_stru
        gen double `u2'=runiform()
        bysort year (`u2'): replace `rank'=_n
        save `draw', replace

        use `main', clear
        merge m:1 year `rank' using `draw', nogen keep(match)
        quietly reghdfe abs_gap placebo_demo placebo_L1_stru $CTRL, ///
            absorb(id province_year) vce(cluster id)
        local bd=_b[placebo_demo]
        local sd=_se[placebo_demo]
        local td=`bd'/`sd'
        local pd=2*ttail(e(df_r),abs(`td'))
        local bs=_b[placebo_L1_stru]
        local ss=_se[placebo_L1_stru]
        local ts=`bs'/`ss'
        local ps=2*ttail(e(df_r),abs(`ts'))
        local joint=(`bd'<0 & `pd'<.10 & `bs'<0 & `ps'<.10)
        post `hyear' (`r') (`bd') (`sd') (`td') (`pd') ///
            (`bs') (`ss') (`ts') (`ps') (`joint')
    restore
    if mod(`r',100)==0 display as text "年份内置换：`r'/1000"
}
postclose `hyear'

use `draws_year', clear
export delimited using "$OUTDIR\年份内随机置换_t值分布_1000次.csv", replace
gen byte abs_demo=abs(t_demo)>=abs(scalar(t_demo_real))
gen byte left_demo=t_demo<=scalar(t_demo_real)
gen byte abs_stru=abs(t_stru)>=abs(scalar(t_stru_real))
gen byte left_stru=t_stru<=scalar(t_stru_real)
quietly count if abs_demo
local pdtwo=(r(N)+1)/($REPS+1)
quietly count if left_demo
local pdleft=(r(N)+1)/($REPS+1)
quietly count if abs_stru
local pstwo=(r(N)+1)/($REPS+1)
quietly count if left_stru
local psleft=(r(N)+1)/($REPS+1)
quietly summarize both_neg10
local bothshare=r(mean)
post `hsum' ("年份内随机置换") ($REPS) ///
    (scalar(b_demo_real)) (scalar(t_demo_real)) (`pdtwo') (`pdleft') ///
    (scalar(b_stru_real)) (scalar(t_stru_real)) (`pstwo') (`psleft') (`bothshare')

histogram t_demo, normal frequency xline(`=scalar(t_demo_real)', ///
    lpattern(dash) lwidth(medthick)) ///
    title("年份内随机置换：demo t值（1000次）") ///
    xtitle("随机置换回归的 demo t 值") ytitle("频数")
graph export "$OUTDIR\图1_年份内随机置换_demo_t值_1000次.png", replace width(2400)
histogram t_stru, normal frequency xline(`=scalar(t_stru_real)', ///
    lpattern(dash) lwidth(medthick)) ///
    title("年份内随机置换：L1_stru t值（1000次）") ///
    xtitle("随机置换回归的 L1_stru t 值") ytitle("频数")
graph export "$OUTDIR\图2_年份内随机置换_L1stru_t值_1000次.png", replace width(2400)

* ---------------------- 方案2：城市内部随机错配年份 ----------------------
tempname hcity
postfile `hcity' int rep ///
    double b_demo se_demo t_demo p_demo b_stru se_stru t_stru p_stru ///
    byte both_neg10 using `draws_city', replace
set seed 20261109
forvalues r=1/$REPS {
    preserve
        use `base_panel', clear
        tempfile main draw
        tempvar u rank u2
        gen double `u'=runiform()
        bysort id (`u'): gen long `rank'=_n
        save `main', replace

        keep id `rank' demo L1_stru
        rename demo placebo_demo
        rename L1_stru placebo_L1_stru
        gen double `u2'=runiform()
        bysort id (`u2'): replace `rank'=_n
        save `draw', replace

        use `main', clear
        merge m:1 id `rank' using `draw', nogen keep(match)
        quietly reghdfe abs_gap placebo_demo placebo_L1_stru $CTRL, ///
            absorb(id province_year) vce(cluster id)
        local bd=_b[placebo_demo]
        local sd=_se[placebo_demo]
        local td=`bd'/`sd'
        local pd=2*ttail(e(df_r),abs(`td'))
        local bs=_b[placebo_L1_stru]
        local ss=_se[placebo_L1_stru]
        local ts=`bs'/`ss'
        local ps=2*ttail(e(df_r),abs(`ts'))
        local joint=(`bd'<0 & `pd'<.10 & `bs'<0 & `ps'<.10)
        post `hcity' (`r') (`bd') (`sd') (`td') (`pd') ///
            (`bs') (`ss') (`ts') (`ps') (`joint')
    restore
    if mod(`r',100)==0 display as text "城市内年份错配：`r'/1000"
}
postclose `hcity'

use `draws_city', clear
export delimited using "$OUTDIR\城市内年份错配_t值分布_1000次.csv", replace
gen byte abs_demo=abs(t_demo)>=abs(scalar(t_demo_real))
gen byte left_demo=t_demo<=scalar(t_demo_real)
gen byte abs_stru=abs(t_stru)>=abs(scalar(t_stru_real))
gen byte left_stru=t_stru<=scalar(t_stru_real)
quietly count if abs_demo
local pdtwo=(r(N)+1)/($REPS+1)
quietly count if left_demo
local pdleft=(r(N)+1)/($REPS+1)
quietly count if abs_stru
local pstwo=(r(N)+1)/($REPS+1)
quietly count if left_stru
local psleft=(r(N)+1)/($REPS+1)
quietly summarize both_neg10
local bothshare=r(mean)
post `hsum' ("城市内年份错配") ($REPS) ///
    (scalar(b_demo_real)) (scalar(t_demo_real)) (`pdtwo') (`pdleft') ///
    (scalar(b_stru_real)) (scalar(t_stru_real)) (`pstwo') (`psleft') (`bothshare')

histogram t_demo, normal frequency xline(`=scalar(t_demo_real)', ///
    lpattern(dash) lwidth(medthick)) ///
    title("城市内年份错配：demo t值（1000次）") ///
    xtitle("随机错配回归的 demo t 值") ytitle("频数")
graph export "$OUTDIR\图3_城市内年份错配_demo_t值_1000次.png", replace width(2400)
histogram t_stru, normal frequency xline(`=scalar(t_stru_real)', ///
    lpattern(dash) lwidth(medthick)) ///
    title("城市内年份错配：L1_stru t值（1000次）") ///
    xtitle("随机错配回归的 L1_stru t 值") ytitle("频数")
graph export "$OUTDIR\图4_城市内年份错配_L1stru_t值_1000次.png", replace width(2400)

postclose `hsum'
use `summary', clear
export delimited using "$OUTDIR\安慰剂_t值分布汇总_1000次.csv", replace
list, noobs abbreviate(24)
display as result "完成：两种安慰剂各1000次，t值分布图及CSV已输出到 $OUTDIR"
log close
