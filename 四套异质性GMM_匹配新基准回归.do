* ============================================================================
* 四套异质性GMM方案：匹配新基准回归的因变量、控制变量和共同样本
* 1 创新中心—外围；2 经济发展水平；3 低人均零售额；4 东/中/西/东北。
* 经济发展水平按2006—2008年实际人均GDP城市均值中位数划分。
*
* 估计期2009—2024；部分动态System-GMM；GMM-style统一1—2期并折叠；
* 控制变量和年份虚拟变量IV-style；表格报告demo、L1_stru和GMM诊断。
* 只输出RTF表，不保存估计数据DTA；不覆盖任何已有文件。
* ============================================================================
version 18.0
clear all
set more off
set maxvar 32767

global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global REGIONMAP "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\省份城市id对照.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\异质性（最后）\最后"
capture mkdir "$OUTDIR"

capture which xtabond2
if _rc {
    display as error "未找到xtabond2。"
    exit 199
}
capture which esttab
if _rc {
    display as error "未找到esttab；请先安装estout。"
    exit 199
}

capture log close _all
log using "$OUTDIR\四套异质性GMM_匹配新基准回归.log", text replace

use "$DATA", clear
foreach v in id year pergdp pergdp_real pergdp_real_growth_pct demo stru ///
    gdp stu mobile retail loan z_net_qty {
    capture confirm variable `v'
    if _rc {
        display as error "估计数据缺少变量 `v'。请核对统一1—2期估计数据。"
        exit 111
    }
}
drop if missing(id) | missing(year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)

* 按新基准回归口径构造被解释变量和控制变量。
foreach v in growth_pp L1_stru pop_proxy pc_stu ln_pc_stu ///
    pc_retail pc_retail_10k ln_mobile pc_loan ln_pc_loan {
    capture drop `v'
}
gen double growth_pp=pergdp_real_growth_pct
gen double L1_stru=L.stru
gen double pop_proxy=gdp*10000/pergdp if gdp>0 & pergdp>0
gen double pc_stu=stu/pop_proxy*10000 if stu>0 & pop_proxy>0
gen double ln_pc_stu=ln(pc_stu) if pc_stu>0
gen double pc_retail=retail/pop_proxy*10000 if retail>0 & pop_proxy>0
gen double pc_retail_10k=pc_retail/10000 if pc_retail<.
gen double ln_mobile=ln(mobile) if mobile>0
gen double pc_loan=loan/pop_proxy*10000 if loan>0 & pop_proxy>0
gen double ln_pc_loan=ln(pc_loan) if pc_loan>0
global CORE demo L1_stru
global CTRL_FINAL ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan

* 与全样本模型对齐：锁定2009—2024年完整的278个城市样本。
egen byte miss_base=rowmiss(growth_pp $CORE $CTRL_FINAL)
gen byte complete_base=inrange(year,2009,2024) & miss_base==0
bysort id: egen int n_complete=total(complete_base)
keep if n_complete==16

egen double mean_growth_t=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-mean_growth_t)
gen double L1_gap=L.abs_gap
gen double L2_gap=L2.abs_gap
foreach v in miss_est estimation retail_pc ///
    base_innovation innovation_rank innovation_group ///
    base_economy low_economy base_retail_pc low_retail_pc ///
    province province_clean region4 {
    capture drop `v'
}
egen byte miss_est=rowmiss(abs_gap L1_gap L2_gap demo L1_stru $CTRL_FINAL)
gen byte estimation=inrange(year,2009,2024) & miss_est==0
quietly count if estimation
if r(N)!=4448 {
    display as error "当前共同估计样本为" r(N) "条，与目标4448条不一致；请检查输入数据和缺失值。"
    exit 459
}

* --- 方案1：创新中心—外围（2006—2008年z_net_qty城市均值，排名前50%） ---
* 注意：按此书面口径重跑后，中心/外围估计与用户指定的旧RTF不一致；
* 不硬编码旧系数，需核实旧表实际使用的分组变量/数据版本后方可声称复现。
bysort id: egen double base_innovation= ///
    mean(cond(inrange(year,2006,2008),z_net_qty,.))
preserve
    keep if year==2009 & estimation
    keep id base_innovation
    drop if missing(id,base_innovation)
    sort id
    quietly count
    local n_innov=r(N)
    egen long innovation_rank=rank(base_innovation), unique
    gen byte innovation_group=innovation_rank>floor(`n_innov'/2)
    keep id innovation_group
    tempfile innov_groups
    save `innov_groups', replace
restore
merge m:1 id using `innov_groups', keep(master match) nogen
quietly count if estimation & missing(innovation_group)
if r(N)>0 {
    display as error "创新中心分组有估计样本未匹配。"
    exit 459
}

* --- 方案2：经济发展水平，按2006—2008年实际人均GDP城市均值中位数分组 ---
bysort id: egen double base_economy= ///
    mean(cond(inrange(year,2006,2008),pergdp_real,.))
quietly summarize base_economy if year==2009 & estimation, detail
local med_economy=r(p50)
if missing(`med_economy') {
    display as error "无法计算实际人均GDP分组中位数。"
    exit 459
}
gen byte low_economy=(base_economy<`med_economy') if base_economy<.

* --- 方案3：低人均零售额（与基准回归的人均零售计算口径一致） ---
gen double retail_pc=pc_retail if pc_retail<.
bysort id: egen double base_retail_pc= ///
    mean(cond(inrange(year,2006,2008),retail_pc,.))
quietly summarize base_retail_pc if year==2009 & estimation, detail
local med_retail=r(p50)
if missing(`med_retail') {
    display as error "无法计算人均零售额分组中位数。"
    exit 459
}
gen byte low_retail_pc=(base_retail_pc<`med_retail') if base_retail_pc<.

* --- 方案4：国家统计局四大区域分类，省份由用户城市—省份ID表映射 ---
merge m:1 id using "$REGIONMAP", keepusing(province) keep(master match) nogen
if missing(province[1]) {
    display as error "省份ID映射未成功并入估计数据。"
    exit 459
}
gen str30 province_clean=ustrregexra(ustrtrim(province), ///
    "省|市|壮族自治区|回族自治区|维吾尔自治区|自治区", "")
gen byte region4=.
replace region4=1 if inlist(province_clean,"北京","天津","河北","上海","江苏") | ///
    inlist(province_clean,"浙江","福建","山东","广东","海南")
replace region4=2 if inlist(province_clean,"山西","安徽","江西","河南","湖北","湖南")
replace region4=3 if inlist(province_clean,"内蒙古","广西","重庆","四川","贵州","云南") | ///
    inlist(province_clean,"西藏","陕西","甘肃","青海","宁夏","新疆")
replace region4=4 if inlist(province_clean,"辽宁","吉林","黑龙江")
label define region4lbl 1 "东部" 2 "中部" 3 "西部" 4 "东北", replace
label values region4 region4lbl
quietly count if estimation & missing(region4)
if r(N)>0 {
    display as error "有" r(N) "条估计观测未能归入东中西东北，请核对省份名称。"
    exit 459
}

* 统一的分组部分动态System-GMM估计程序（不加入组别交互项）。
capture program drop fit_subgroup
program define fit_subgroup, eclass
    syntax [if], NAME(name)
    marksample touse
    quietly xtabond2 abs_gap L1_gap L2_gap demo L1_stru $CTRL_FINAL i.year ///
        if `touse' & estimation, ///
        gmmstyle(L1_gap L2_gap demo L1_stru, laglimits(1 2) collapse) ///
        ivstyle($CTRL_FINAL i.year) ///
        twostep robust small orthogonal
    quietly test demo
    estadd scalar Demo_p=r(p)
    quietly test L1_stru
    estadd scalar Stru_p=r(p)
    quietly test L1_gap L2_gap
    estadd scalar DynJoint_p=r(p)
    estadd scalar AR1_p=e(ar1p)
    estadd scalar AR2_p=e(ar2p)
    estadd scalar Hansen_p=e(hansenp)
    estadd scalar Instruments=e(j)
    estadd scalar Groups=e(N_g)
    estimates store `name'
end

local keepvars demo L1_stru
local coeflabs demo "demo" L1_stru "L1_stru"
local statlist N Groups Instruments Demo_p Stru_p DynJoint_p AR1_p AR2_p Hansen_p
local statlabels "观测值" "城市数" "工具变量数" "demo p值" "L1_stru p值" ///
    "动态项联合检验p值" "AR(1) p值" "AR(2) p值" "Hansen p值"

* 方案1：创新中心—外围。
eststo clear
fit_subgroup, name(s1_full)
fit_subgroup if innovation_group==1, name(s1_center)
fit_subgroup if innovation_group==0, name(s1_periphery)
esttab s1_full s1_center s1_periphery ///
    using "$OUTDIR\方案1_创新中心外围.rtf", replace rtf ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) keep(`keepvars') order(`keepvars') ///
    coeflabels(`coeflabs') mtitles("全样本" "创新中心" "创新外围") ///
    stats(`statlist', labels(`statlabels') fmt(0 0 0 4 4 4 4 4 4)) ///
    title("异质性部分动态System-GMM：创新中心—外围") ///
    addnotes("创新中心为2006—2008年z_net_qty城市均值排名前50%的城市。", ///
      "第一列为全样本，后两列为创新中心、创新外围；只展示demo和L1_stru。", ///
      "因变量、四项控制变量和共同样本均与新基准回归一致；GMM设定保持不变。")

* 方案2：经济发展水平（现有对应结果表使用实际人均GDP作分组指标）。
eststo clear
fit_subgroup, name(s2_full)
fit_subgroup if low_economy==1, name(s2_low)
fit_subgroup if low_economy==0, name(s2_high)
esttab s2_full s2_low s2_high ///
    using "$OUTDIR\方案2_经济发展水平高低.rtf", replace rtf ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) keep(`keepvars') order(`keepvars') ///
    coeflabels(`coeflabs') mtitles("全样本" "低经济发展组" "高经济发展组") ///
    stats(`statlist', labels(`statlabels') fmt(0 0 0 4 4 4 4 4 4)) ///
    title("异质性部分动态System-GMM：经济发展水平高低") ///
    addnotes("经济发展水平按2006—2008年实际人均GDP城市均值中位数划分。", ///
      "第一列为全样本；表格只展示demo和L1_stru。", ///
      "因变量、四项控制变量和共同样本均与新基准回归一致；GMM设定保持不变。")

* 方案3：低人均零售额组（零售额/估算人口×10000）。
eststo clear
fit_subgroup, name(s3_full)
fit_subgroup if low_retail_pc==1, name(s3_low)
fit_subgroup if low_retail_pc==0, name(s3_high)
esttab s3_full s3_low s3_high ///
    using "$OUTDIR\方案3_低人均零售低高组.rtf", replace rtf ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) keep(`keepvars') order(`keepvars') ///
    coeflabels(`coeflabs') mtitles("全样本" "低人均零售组" "高人均零售组") ///
    stats(`statlist', labels(`statlabels') fmt(0 0 0 4 4 4 4 4 4)) ///
    title("异质性部分动态System-GMM：人均零售额") ///
    addnotes("人均零售额=社会消费品零售总额/估算人口×10000；按2006—2008年城市均值中位数划分。", ///
      "分组指标不是零售缺口，也不是零售额/GDP。", ///
      "因变量、四项控制变量和共同样本均与新基准回归一致；GMM设定保持不变。")

* 方案4：东、中、西、东北四大区域；采用国家统计局地区划分。
eststo clear
fit_subgroup, name(s4_full)
fit_subgroup if region4==1, name(s4_east)
fit_subgroup if region4==2, name(s4_central)
fit_subgroup if region4==3, name(s4_west)
fit_subgroup if region4==4, name(s4_northeast)
esttab s4_full s4_east s4_central s4_west s4_northeast ///
    using "$OUTDIR\方案4_东中西东北四大区域.rtf", replace rtf ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) keep(`keepvars') order(`keepvars') ///
    coeflabels(`coeflabs') mtitles("全样本" "东部" "中部" "西部" "东北") ///
    stats(`statlist', labels(`statlabels') fmt(0 0 0 4 4 4 4 4 4)) ///
    title("异质性部分动态System-GMM：东中西东北地区") ///
    addnotes("按国家统计局东、中、西、东北四大区域口径划分。", ///
      "第一列为全样本；后四列依次为东部、中部、西部和东北。", ///
      "因变量、四项控制变量和共同样本均与新基准回归一致；GMM设定保持不变。")

display as result "四套异质性方案估计和RTF表格输出完成；未保存任何DTA文件。"
log close
