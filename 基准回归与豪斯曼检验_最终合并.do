version 18.0
clear all
set more off
set maxvar 32767

* 区域经济差异：阶梯式基准回归、Hausman检验及统计量小表
global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global REGIONMAP "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\省份城市id对照.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\基准回归（最后）\最后"
capture mkdir "$OUTDIR"
capture log close _all
log using "$OUTDIR\基准回归与豪斯曼检验_最终合并.log", text replace

foreach cmd in reghdfe esttab eststo estadd {
    capture which `cmd'
    if _rc exit 199
}

* 一、统一整理数据和变量
use "$DATA", clear
drop if missing(id) | missing(year)
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
gen double pc_stu=stu/pop_proxy*10000 if stu>0 & pop_proxy>0
gen double ln_pc_stu=ln(pc_stu) if pc_stu>0
gen double pc_retail=retail/pop_proxy*10000 if retail>0 & pop_proxy>0
gen double pc_retail_10k=pc_retail/10000 if pc_retail<.
gen double ln_mobile=ln(mobile) if mobile>0
gen double pc_loan=loan/pop_proxy*10000 if loan>0 & pop_proxy>0
gen double ln_pc_loan=ln(pc_loan) if pc_loan>0

merge m:1 id using "$REGIONMAP", keepusing(province) keep(master match) nogen
encode province, gen(province_id)
egen long province_year=group(province_id year)

egen byte miss_final=rowmiss(abs_gap demo L1_stru ln_pc_stu ///
    pc_retail_10k ln_mobile ln_pc_loan province_id)
gen byte valid_final=inrange(year,2009,2024) & miss_final==0
bysort id: egen int n_valid_final=total(valid_final)
keep if n_valid_final==16
keep if valid_final==1
xtset id year

global X demo L1_stru ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan

* 二、阶梯式双固定效应基准回归
eststo clear
eststo M1: reghdfe abs_gap demo L1_stru, ///
    absorb(id province_year) vce(cluster id)
eststo M2: reghdfe abs_gap demo L1_stru ln_pc_stu, ///
    absorb(id province_year) vce(cluster id)
eststo M3: reghdfe abs_gap demo L1_stru ln_pc_stu pc_retail_10k, ///
    absorb(id province_year) vce(cluster id)
eststo M4: reghdfe abs_gap demo L1_stru ln_pc_stu pc_retail_10k ln_mobile, ///
    absorb(id province_year) vce(cluster id)
eststo M5: reghdfe abs_gap $X, ///
    absorb(id province_year) vce(cluster id)

esttab M1 M2 M3 M4 M5 using "$OUTDIR\表1_阶梯式基准回归.rtf", ///
    replace rtf b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(demo L1_stru ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan) ///
    order(demo L1_stru ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan) ///
    mtitles("(1) 核心变量" "(2) +人力资本" "(3) +零售规模" ///
        "(4) +移动电话" "(5) +人均贷款") ///
    stats(N r2, fmt(0 4) labels("观测值" "总体 R-squared (reghdfe)")) ///
    title("区域经济差异：reghdfe 阶梯式双固定效应基准回归") ///
    addnotes("因变量为城市实际人均GDP增长率相对当年城市均值的绝对偏离；数值越小表示差异越小。", ///
        "人力资本为人均普通高校在校生数的对数；人均零售按万元/人计量；移动电话为年末用户总量的对数；人均贷款取对数。", ///
        "数据为2009—2024年278个城市的4448个观测值，所有列样本一致。", ///
        "控制城市固定效应和省份×年份固定效应，标准误按城市聚类。")

* 三、经典Hausman检验
eststo clear
quietly xtreg abs_gap $X i.province_year, fe
estimates store FE_classic
quietly xtreg abs_gap $X i.province_year, re
estimates store RE_classic
hausman FE_classic RE_classic, sigmamore
local H_chi2=r(chi2)
local H_df=r(df)
local H_p=r(p)

* 四、城市聚类稳健的Mundlak-Hausman检验
bysort id: egen double mean_demo=mean(demo)
bysort id: egen double mean_L1_stru=mean(L1_stru)
bysort id: egen double mean_ln_pc_stu=mean(ln_pc_stu)
bysort id: egen double mean_pc_retail_10k=mean(pc_retail_10k)
bysort id: egen double mean_ln_mobile=mean(ln_mobile)
bysort id: egen double mean_ln_pc_loan=mean(ln_pc_loan)

quietly xtreg abs_gap $X mean_demo mean_L1_stru mean_ln_pc_stu ///
    mean_pc_retail_10k mean_ln_mobile mean_ln_pc_loan ///
    i.province_year, re vce(cluster id)
test mean_demo mean_L1_stru mean_ln_pc_stu mean_pc_retail_10k ///
    mean_ln_mobile mean_ln_pc_loan
local MH_chi2=r(chi2)
local MH_df=r(df)
local MH_p=r(p)

* 五、详细FE与RE比较表
quietly xtreg abs_gap $X i.province_year, fe vce(cluster id)
estadd scalar Classical_chi2=`H_chi2'
estadd scalar Classical_p=`H_p'
estadd scalar Mundlak_chi2=`MH_chi2'
estadd scalar Mundlak_p=`MH_p'
estimates store FE_robust

quietly xtreg abs_gap $X i.province_year, re vce(cluster id)
estadd scalar Classical_chi2=`H_chi2'
estadd scalar Classical_p=`H_p'
estadd scalar Mundlak_chi2=`MH_chi2'
estadd scalar Mundlak_p=`MH_p'
estimates store RE_robust

esttab FE_robust RE_robust using "$OUTDIR\表2_固定效应与随机效应比较.rtf", ///
    replace rtf b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(demo L1_stru ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan) ///
    order(demo L1_stru ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan) ///
    mtitles("固定效应 FE" "随机效应 RE") ///
    stats(N Classical_chi2 Classical_p Mundlak_chi2 Mundlak_p, ///
        fmt(0 4 4 4 4) labels("观测值" "经典Hausman chi2" ///
        "经典Hausman p值" "稳健Mundlak chi2" "稳健Mundlak p值")) ///
    title("区域经济差异：固定效应与随机效应比较") ///
    addnotes("P值小于0.05时拒绝随机效应，选择固定效应模型。", ///
        "经典Hausman因协方差差矩阵非正定仅作参考，主要依据城市聚类稳健的Mundlak-Hausman检验。")

* 六、按论文常用样式输出单行Hausman统计量表
matrix HAUSMAN_MAIN=(`MH_chi2',`MH_df',`MH_p')
matrix rownames HAUSMAN_MAIN=Hausman
matrix colnames HAUSMAN_MAIN=Chi2 df p_value

esttab matrix(HAUSMAN_MAIN, fmt(4 0 4)) using ///
    "$OUTDIR\表3_Hausman统计量汇总.rtf", ///
    replace rtf nomtitles nonumber ///
    collabels("Chi2" "自由度" "P值") ///
    title("Hausman 统计量：最终控制变量") ///
    addnotes("本表报告城市聚类稳健的Mundlak-Hausman检验。", ///
        "原假设为随机效应模型一致；P值小于0.05时选择固定效应模型。")

display as result "Classical Hausman chi2(" `H_df' ") = " %9.4f `H_chi2'
display as result "Classical Hausman p-value = " %9.6f `H_p'
display as result "Robust Mundlak chi2(" `MH_df' ") = " %9.4f `MH_chi2'
display as result "Robust Mundlak p-value = " %9.6f `MH_p'

log close
