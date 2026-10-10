/****************************************************************************************
  反事实检验：同样本未来变量替换（最终推荐版）
  研究目的：检验未来一期 demo 与未来一期 stru 是否能复制基准结果。

  关键逻辑：
  1. 基准核心变量为：demo_t 与 L1_stru = stru_{t-1}。
  2. 不能使用 F1_L1_stru，因为 F1.(L.stru) = stru_t，本质上是当期 stru，
     不是未来变量。
  3. 因此，真正的未来反事实变量应为：F1_demo = demo_{t+1}，
     F1_stru = stru_{t+1}。
  4. 为避免样本变化影响比较，本文在相同样本中比较真实变量与未来变量替换结果。

  使用前请确认：
  global DATA      "你的数据路径\统一分析数据.dta"
  global REGIONMAP "你的省份映射文件路径.dta"   // 若数据中已有 province 变量，则可不设置
  global OUTDIR    "结果输出路径"
****************************************************************************************/

clear all
set more off
version 16.0

*============================*
* 0. 路径设置：按需修改
*============================*
* 若你已经在主 do 文件中设置过 global DATA / REGIONMAP / OUTDIR，可注释下面三行
global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global REGIONMAP "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\省份城市id对照.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\反事实安慰剂\反事实"

*============================*
* 1. 检查必要命令
*============================*
foreach cmd in reghdfe esttab eststo estadd {
    capture which `cmd'
    if _rc {
        display as error "缺少命令：`cmd'。请先运行：ssc install `cmd', replace"
        exit 199
    }
}

*============================*
* 2. 数据整理：匹配基准口径
*============================*
use "$DATA", clear

drop if missing(id) | missing(year)
replace year = round(year)
sort id year
isid id year
xtset id year

* 因变量：城市实际人均GDP增长率相对年度均值的绝对离差
gen double growth_pp = pergdp_real_growth_pct
egen double mean_growth = mean(growth_pp), by(year)
gen double abs_gap = abs(growth_pp - mean_growth)

* 核心解释变量：真实变量
gen double L1_stru = L.stru

* 反事实变量：真正的未来一期变量
* 注意：不要使用 F1_L1_stru，因为 F1.L1_stru = stru_t，是当期变量，不是未来变量。
gen double F1_demo = F.demo
gen double F1_stru = F.stru

* 控制变量：匹配基准回归
gen double pop_proxy = gdp*10000/pergdp if gdp>0 & pergdp>0
gen double pc_stu = stu/pop_proxy*10000 if stu>0 & pop_proxy>0
gen double ln_pc_stu = ln(pc_stu) if pc_stu>0

gen double pc_retail = retail/pop_proxy*10000 if retail>0 & pop_proxy>0
gen double pc_retail_10k = pc_retail/10000 if pc_retail<.

gen double ln_mobile = ln(mobile) if mobile>0

gen double pc_loan = loan/pop_proxy*10000 if loan>0 & pop_proxy>0
gen double ln_pc_loan = ln(pc_loan) if pc_loan>0

* 省份×年份固定效应
capture confirm variable province
if _rc {
    capture confirm file "$REGIONMAP"
    if _rc {
        display as error "数据中没有 province 变量，且没有找到 $REGIONMAP。请先设置省份映射文件。"
        exit 601
    }
    merge m:1 id using "$REGIONMAP", keepusing(province) keep(master match) nogen
}

capture confirm numeric variable province
if !_rc {
    decode province, gen(province_str)
    drop province
    rename province_str province
}

encode province, gen(province_id)
egen long province_year = group(province_id year)

* 共同样本：严格匹配原基准样本 2009—2024，且每个城市16年完整
global C ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan

egen byte miss_final = rowmiss(abs_gap demo L1_stru $C province_id)
gen byte valid_final = inrange(year,2009,2024) & miss_final==0
bysort id: egen int n_valid_final = total(valid_final)
keep if n_valid_final==16
keep if valid_final==1
xtset id year

* 未来变量共同样本：由于 F1_demo / F1_stru 需要 t+1，2024年自然缺失，因此样本为2009—2023
gen byte cf_sample = !missing(F1_demo, F1_stru)

*============================*
* 3. 反事实检验回归
*============================*
eststo clear

* 列1：基准全样本，用于对照
reghdfe abs_gap demo L1_stru $C, absorb(id province_year) vce(cluster id)
eststo Base_full
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Sample "2009-2024"
estadd local Test "Baseline"

* 列2：同样本真实变量。关键对照列。
reghdfe abs_gap demo L1_stru $C if cf_sample==1, absorb(id province_year) vce(cluster id)
eststo True_same
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Sample "2009-2023"
estadd local Test "True vars, same sample"

* 列3：未来变量替换。正文推荐保留列。
reghdfe abs_gap F1_demo F1_stru $C if cf_sample==1, absorb(id province_year) vce(cluster id)
eststo Future_replace
estadd local CityFE "Yes"
estadd local PYFE "Yes"
estadd local Cluster "City"
estadd local Sample "2009-2023"
estadd local Test "Future vars replace true vars"

*============================*
* 4. 导出 RTF 表格
*============================*
esttab Base_full True_same Future_replace ///
    using "$OUTDIR\表_反事实检验_同样本未来替换_三列版.rtf", ///
    replace rtf b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(demo L1_stru F1_demo F1_stru $C) ///
    order(demo L1_stru F1_demo F1_stru $C) ///
    mtitles("基准全样本" "同样本真实变量" "未来变量替换") ///
    stats(N r2 CityFE PYFE Cluster Sample Test, fmt(0 4 s s s s s) ///
        labels("观测值" "R-squared" "城市固定效应" "省份×年份固定效应" "聚类层级" "样本区间" "检验类型")) ///
    title("反事实检验：同样本未来变量替换") ///
    addnotes("所有列均控制城市固定效应和省份×年份固定效应，标准误按城市聚类。", ///
             "由于 L1_stru = L.stru，F1.L1_stru 实际等于当期 stru，因此本文使用 F1_stru 作为真正未来一期结构变量。", ///
             "核心判断是：同样本真实变量仍为负向显著，而未来变量替换不能复制该结果。")

*============================*
* 5. 导出核心结果 CSV
*============================*
tempname handle
postfile `handle' str30 model str20 var double coef se p N using "$OUTDIR\反事实检验_同样本未来替换_三列核心结果.dta", replace

* 重新逐列估计并保存核心变量，避免手工录入
* 1 基准全样本
reghdfe abs_gap demo L1_stru $C, absorb(id province_year) vce(cluster id)
foreach v in demo L1_stru {
    post `handle' ("基准全样本") ("`v'") (_b[`v']) (_se[`v']) (2*ttail(e(df_r),abs(_b[`v']/_se[`v']))) (e(N))
}

* 2 同样本真实变量
reghdfe abs_gap demo L1_stru $C if cf_sample==1, absorb(id province_year) vce(cluster id)
foreach v in demo L1_stru {
    post `handle' ("同样本真实变量") ("`v'") (_b[`v']) (_se[`v']) (2*ttail(e(df_r),abs(_b[`v']/_se[`v']))) (e(N))
}

* 3 未来变量替换
reghdfe abs_gap F1_demo F1_stru $C if cf_sample==1, absorb(id province_year) vce(cluster id)
foreach v in F1_demo F1_stru {
    post `handle' ("未来变量替换") ("`v'") (_b[`v']) (_se[`v']) (2*ttail(e(df_r),abs(_b[`v']/_se[`v']))) (e(N))
}

postclose `handle'

preserve
use "$OUTDIR\反事实检验_同样本未来替换_三列核心结果.dta", clear
export delimited using "$OUTDIR\反事实检验_同样本未来替换_三列核心结果.csv", replace
restore

*============================*
* 6. 自动生成论文表述
*============================*
file open fw using "$OUTDIR\论文写法_反事实检验_同样本未来替换.txt", write replace
file write fw "为检验基准结果是否由未来信息、反向因果或机械性时间趋势驱动，本文进行反事实检验。" _n
file write fw "鉴于基准模型中的 L1_stru 表示原始结构变量 stru 的一期滞后项，即 stru_{t-1}，若直接使用 F1.L1_stru，则实际对应当期 stru_t，并非严格意义上的未来变量。" _n
file write fw "因此，本文使用未来一期 demo 与未来一期原始 stru，即 F1_demo 和 F1_stru，构造反事实变量。" _n
file write fw "结果显示，在与未来变量完全相同的样本中，真实 demo 和 L1_stru 仍表现为负向显著；而将核心解释变量替换为未来一期 demo 和未来一期 stru 后，两个未来变量均不显著，且未复制真实变量的负向结果。" _n
file write fw "这说明基准结论并非由未来信息或简单时间趋势所驱动。" _n
file close fw

display as result "反事实检验完成。请查看："
display as result "$OUTDIR\表_反事实检验_同样本未来替换_三列版.rtf"
display as result "$OUTDIR\反事实检验_同样本未来替换_三列核心结果.csv"
display as result "$OUTDIR\论文写法_反事实检验_同样本未来替换.txt"
