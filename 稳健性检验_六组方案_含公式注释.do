/****************************************************************************************
  稳健性检验：六组方案（3个替代被解释变量 + 剔除2020 + 两种5%缩尾）

  目标：在同一张表中保留六组稳健性检验：
  (1) 年度秩次绝对离差 ×100
  (2) 实际增长率标准化绝对离差 ×100
  (3) 年度正态秩转换绝对离差 ×100
  (4) 剔除2020年
  (5) abs_gap全样本5%缩尾 ×100
  (6) abs_gap按年份5%缩尾 ×100

  说明：
  1. 前三列为替换被解释变量；后三列为样本/极端值处理稳健性。
  2. 所有列均控制城市固定效应和省份×年份固定效应，标准误按城市聚类。
  3. 因变量乘以100仅改变计量单位，不改变t值、p值和显著性。
****************************************************************************************/

clear all
set more off

*===========================*
* 一、路径设置：运行前按需修改
*===========================*
global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global REGIONMAP "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\省份城市id对照.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\稳健性检验"

capture mkdir "$OUTDIR"

*===========================*
* 二、依赖命令检查
*===========================*
foreach cmd in reghdfe esttab eststo estadd {
    capture which `cmd'
    if _rc {
        display as error "缺少命令：`cmd'。请先运行：ssc install `cmd', replace"
        exit 199
    }
}

*===========================*
* 三、数据整理：严格匹配基准样本
*===========================*
use "$DATA", clear

drop if missing(id) | missing(year)
replace year = round(year)
sort id year
isid id year
xtset id year

* 基准增长率与基准因变量
* growth_pp：实际人均GDP增长率，单位为百分点
* abs_gap_it = | growth_pp_it - 年度平均growth_pp_t |
gen double growth_pp = pergdp_real_growth_pct
egen double mean_growth = mean(growth_pp), by(year)
gen double abs_gap = abs(growth_pp - mean_growth)
gen double abs_gap100 = abs_gap * 100
label var abs_gap "年度均值绝对离差"
label var abs_gap100 "年度均值绝对离差×100"

* 核心解释变量
gen double L1_stru = L.stru

* 基准控制变量
gen double pop_proxy = gdp*10000/pergdp if gdp>0 & pergdp>0
gen double pc_stu = stu/pop_proxy*10000 if stu>0 & pop_proxy>0
gen double ln_pc_stu = ln(pc_stu) if pc_stu>0
gen double pc_retail = retail/pop_proxy*10000 if retail>0 & pop_proxy>0
gen double pc_retail_10k = pc_retail/10000 if pc_retail<.
gen double ln_mobile_base = ln(mobile) if mobile>0
gen double pc_loan = loan/pop_proxy*10000 if loan>0 & pop_proxy>0
gen double ln_pc_loan = ln(pc_loan) if pc_loan>0

global C ln_pc_stu pc_retail_10k ln_mobile_base ln_pc_loan

* 省份×年份固定效应
capture confirm variable province
if _rc {
    merge m:1 id using "$REGIONMAP", keepusing(province) keep(master match) nogen
}
encode province, gen(province_id)
egen long province_year = group(province_id year)

* 共同样本：2009—2024，每个城市16年完整，匹配原基准口径
egen byte miss_final = rowmiss(abs_gap demo L1_stru $C province_id)
gen byte valid_final = inrange(year,2009,2024) & miss_final==0
bysort id: egen int n_valid_final = total(valid_final)
keep if n_valid_final==16
keep if valid_final==1
sort id year
xtset id year

*===========================*
* 四、构造三组替代被解释变量
*===========================*

* 方案1：年度秩次绝对离差 ×100
* 公式：Y1_it = | ((rank(growth_pp_it)-0.5)/N_t) - 0.5 | × 100
* 含义：城市增长率在当年城市分布中的相对位置，距离中位秩次越远，差异越大。
bysort year: egen double rank_growth = rank(growth_pp)
bysort year: egen double n_growth = count(growth_pp)
gen double rankpct_growth = (rank_growth - 0.5) / n_growth
gen double y_rank_abs100 = abs(rankpct_growth - 0.5) * 100
label var y_rank_abs100 "年度秩次绝对离差×100"

* 方案2：实际增长率标准化绝对离差 ×100
* 公式：Y2_it = | (growth_pp_it - mean_t(growth_pp)) / sd_t(growth_pp) | × 100
* 含义：按年度均值和年度标准差标准化后，度量城市增长率偏离年度中心水平的程度。
egen double mean_growth_pp = mean(growth_pp), by(year)
egen double sd_growth_pp = sd(growth_pp), by(year)
gen double y_std_growth_abs100 = abs((growth_pp - mean_growth_pp) / sd_growth_pp) * 100 if sd_growth_pp>0
label var y_std_growth_abs100 "实际增长率标准化绝对离差×100"

* 方案3：年度正态秩转换绝对离差 ×100
* 公式：Y3_it = | Φ^(-1)( (rank(growth_pp_it)-0.5)/N_t ) | × 100
* 含义：先把年度增长率秩次转换为近似正态分位数，再取绝对离差，对极端值和原始量纲更不敏感。
gen double p_rank_growth = rankpct_growth
replace p_rank_growth = 0.0001 if p_rank_growth<=0
replace p_rank_growth = 0.9999 if p_rank_growth>=1
gen double y_normrank_abs100 = abs(invnormal(p_rank_growth)) * 100
label var y_normrank_abs100 "年度正态秩转换绝对离差×100"

*===========================*
* 五、构造两种5%缩尾因变量
*===========================*

* 方案5：abs_gap全样本5%缩尾 ×100
* 设定：用全样本abs_gap的5%和95%分位数作为上下界进行缩尾。
quietly _pctile abs_gap, p(5 95)
scalar p5_all = r(r1)
scalar p95_all = r(r2)
gen double abs_gap_w5_all = abs_gap
replace abs_gap_w5_all = p5_all if abs_gap_w5_all < p5_all & abs_gap_w5_all < .
replace abs_gap_w5_all = p95_all if abs_gap_w5_all > p95_all & abs_gap_w5_all < .
gen double abs_gap_w5_all100 = abs_gap_w5_all * 100
label var abs_gap_w5_all100 "全样本5%缩尾绝对离差×100"

* 方案6：abs_gap按年份5%缩尾 ×100
* 设定：每一年单独计算abs_gap的5%和95%分位数，并在年度内部进行缩尾。
gen double abs_gap_w5_year = abs_gap
levelsof year, local(years)
foreach yy of local years {
    quietly _pctile abs_gap if year == `yy', p(5 95)
    local p5_y = r(r1)
    local p95_y = r(r2)
    replace abs_gap_w5_year = `p5_y' if year == `yy' & abs_gap_w5_year < `p5_y' & abs_gap_w5_year < .
    replace abs_gap_w5_year = `p95_y' if year == `yy' & abs_gap_w5_year > `p95_y' & abs_gap_w5_year < .
}
gen double abs_gap_w5_year100 = abs_gap_w5_year * 100
label var abs_gap_w5_year100 "按年份5%缩尾绝对离差×100"

*===========================*
* 六、六组稳健性回归
*===========================*
eststo clear

* (1) 替代被解释变量：年度秩次绝对离差×100
reghdfe y_rank_abs100 demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo R1_rank
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "替换被解释变量"
estadd local Setting "秩次离差"
estadd local Unit "×100"

* (2) 替代被解释变量：实际增长率标准化绝对离差×100
reghdfe y_std_growth_abs100 demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo R2_stdgrowth
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "替换被解释变量"
estadd local Setting "标准化离差"
estadd local Unit "×100"

* (3) 替代被解释变量：年度正态秩转换绝对离差×100
reghdfe y_normrank_abs100 demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo R3_normrank
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "替换被解释变量"
estadd local Setting "正态秩离差"
estadd local Unit "×100"

* (4) 剔除2020年：使用基准因变量abs_gap×100
reghdfe abs_gap100 demo L1_stru $C if year != 2020, absorb(id province_year) vce(cluster id)
eststo R4_drop2020
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "样本调整"
estadd local Setting "剔除2020年"
estadd local Unit "×100"

* (5) 全样本5%缩尾：使用abs_gap全样本5%/95%缩尾后×100
reghdfe abs_gap_w5_all100 demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo R5_w5all
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "极端值处理"
estadd local Setting "全样本5%缩尾"
estadd local Unit "×100"

* (6) 按年份5%缩尾：使用abs_gap年度内部5%/95%缩尾后×100
reghdfe abs_gap_w5_year100 demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo R6_w5year
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Type "极端值处理"
estadd local Setting "按年份5%缩尾"
estadd local Unit "×100"

*===========================*
* 七、导出RTF表格：表格下方附每个方案的设定与公式
*===========================*
esttab R1_rank R2_stdgrowth R3_normrank R4_drop2020 R5_w5all R6_w5year ///
    using "$OUTDIR\表_稳健性检验_六组方案.rtf", ///
    replace rtf b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(demo L1_stru $C) ///
    order(demo L1_stru $C) ///
    mtitles("秩次离差×100" "实际增长标准化离差×100" "正态秩离差×100" "剔除2020年" "全样本5%缩尾" "按年份5%缩尾") ///
    stats(N r2 CityFE PYFE Cluster Type Setting Unit, fmt(0 4 s s s s s s) ///
        labels("观测值" "R-squared" "城市固定效应" "省份×年份固定效应" "聚类层级" "检验类型" "方案设定" "因变量单位")) ///
    title("稳健性检验：六组方案") ///
    addnotes("所有列均控制城市固定效应和省份×年份固定效应，标准误按城市聚类。", ///
             "方案1 秩次离差：Y=|((rank(growth_pp)-0.5)/N_t)-0.5|×100，衡量城市增长率在年度分布中偏离中位秩次的程度。", ///
             "方案2 实际增长标准化离差：Y=|(growth_pp-mean_t(growth_pp))/sd_t(growth_pp)|×100，按年度波动幅度标准化增长率偏离。", ///
             "方案3 正态秩离差：Y=|invnormal((rank(growth_pp)-0.5)/N_t)|×100，将年度秩次转换为正态分位数后度量离差。", ///
             "方案4 剔除2020年：使用基准因变量abs_gap×100，但删除2020年样本，以排除特殊年份冲击。", ///
             "方案5 全样本5%缩尾：对全样本abs_gap按5%和95%分位数缩尾后乘以100。", ///
             "方案6 按年份5%缩尾：在每一年内部对abs_gap按5%和95%分位数缩尾后乘以100。", ///
             "注：growth_pp为实际人均GDP增长率；abs_gap=|growth_pp-年度均值|；因变量乘以100只改变计量单位，不改变t值和显著性。")

*===========================*
* 八、导出核心结果CSV
*===========================*
tempname handle
postfile `handle' str40 scheme str20 var double coef se p N using "$OUTDIR\稳健性检验_六组方案_核心结果.dta", replace

foreach m in R1_rank R2_stdgrowth R3_normrank R4_drop2020 R5_w5all R6_w5year {
    est restore `m'
    local n = e(N)
    foreach v in demo L1_stru {
        local b = _b[`v']
        local s = _se[`v']
        local p = 2*ttail(e(df_r), abs(_b[`v']/_se[`v']))
        post `handle' ("`m'") ("`v'") (`b') (`s') (`p') (`n')
    }
}
postclose `handle'

preserve
use "$OUTDIR\稳健性检验_六组方案_核心结果.dta", clear
export delimited using "$OUTDIR\稳健性检验_六组方案_核心结果.csv", replace
restore

*===========================*
* 九、另存一份方案说明TXT，方便复制到论文正文或表下注释
*===========================*
file open note using "$OUTDIR\稳健性检验_六组方案_表下注释.txt", write replace
file write note "稳健性检验六组方案说明" _n
file write note "方案1：秩次离差×100。公式：Y=|((rank(growth_pp)-0.5)/N_t)-0.5|×100；衡量城市增长率在当年城市分布中偏离中位秩次的程度。" _n
file write note "方案2：实际增长标准化离差×100。公式：Y=|(growth_pp-mean_t(growth_pp))/sd_t(growth_pp)|×100；按年度标准差标准化城市增长率偏离。" _n
file write note "方案3：正态秩离差×100。公式：Y=|invnormal((rank(growth_pp)-0.5)/N_t)|×100；将年度秩次转换为正态分位数后取绝对离差。" _n
file write note "方案4：剔除2020年。使用基准因变量abs_gap×100，删除2020年样本，以排除特殊年份冲击。" _n
file write note "方案5：全样本5%缩尾。对全样本abs_gap按5%和95%分位数缩尾后乘以100。" _n
file write note "方案6：按年份5%缩尾。在每一年内部对abs_gap按5%和95%分位数缩尾后乘以100。" _n
file write note "其中 growth_pp 为实际人均GDP增长率，abs_gap=|growth_pp-年度均值|。" _n
file close note

display as text "完成：已输出 表_稳健性检验_六组方案.rtf、核心结果CSV和表下注释TXT。"
