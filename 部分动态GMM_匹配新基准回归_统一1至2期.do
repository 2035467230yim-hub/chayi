* 区域经济差异：部分动态System-GMM（匹配新基准回归控制变量）
* 所有GMM-style变量统一使用1—2期滞后工具并折叠。
version 18.0
clear all
set more off
set maxvar 32767
global DATA "D:\Desktop\工作\研究生师门工作\专利转移\【最终成功】专利转移与区域经济增长【2009-2024年】\删除部分城市和缩尾1%的数据.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（增长率绝对离差）\GMM（最后）\最后"
global OUTLOG "$OUTDIR\部分动态GMM_匹配新基准回归_统一1至2期.log"
global OUTTAB "$OUTDIR\表_部分动态GMM_匹配新基准回归_统一1至2期.rtf"
capture mkdir "$OUTDIR"
capture log close _all
log using "$OUTLOG", text replace
capture which xtabond2
if _rc {
    display as error "未找到xtabond2"
    exit 199
}
capture which esttab
if _rc {
    display as error "未找到esttab"
    exit 199
}
use "$DATA", clear
drop if missing(id) | missing(year)
replace year=round(year)
sort id year
isid id year
xtset id year
keep if inrange(year,2005,2024)
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
egen byte miss_base=rowmiss(growth_pp $CORE $CTRL_FINAL)
gen byte complete_base=inrange(year,2009,2024) & miss_base==0
bysort id: egen int n_complete=total(complete_base)
keep if n_complete==16
egen double mean_growth_t=mean(growth_pp), by(year)
gen double abs_gap=abs(growth_pp-mean_growth_t)
gen double L1_gap=L.abs_gap
gen double L2_gap=L2.abs_gap
egen byte miss_est=rowmiss(abs_gap L1_gap L2_gap $CORE $CTRL_FINAL)
gen byte estimation=inrange(year,2009,2024) & miss_est==0
quietly count if estimation
local N_est=r(N)
egen byte citytag=tag(id) if estimation
quietly count if citytag
local N_city=r(N)
display as result "统一1—2期：估计观测值 N=`N_est'，城市数 G=`N_city'"
xtabond2 abs_gap L1_gap L2_gap $CORE $CTRL_FINAL i.year if estimation, ///
    gmmstyle(L1_gap L2_gap demo L1_stru, laglimits(1 2) collapse) ///
    ivstyle($CTRL_FINAL i.year) ///
    twostep robust small orthogonal
eststo MAIN
test L1_gap L2_gap
estadd scalar DynJoint_p=r(p)
lincom L1_gap+L2_gap
estadd scalar DynSum=r(estimate)
estadd scalar AR1_p=e(ar1p)
estadd scalar AR2_p=e(ar2p)
estadd scalar Hansen_p=e(hansenp)
estadd scalar Instruments=e(j)
estadd scalar Groups=e(N_g)
display as result "动态项联合p值=" e(DynJoint_p)
display as result "动态系数和=" e(DynSum)
display as result "AR(1) p值=" e(ar1p)
display as result "AR(2) p值=" e(ar2p)
display as result "Hansen p值=" e(hansenp)
display as result "工具变量数=" e(j)
esttab MAIN using "$OUTTAB", replace rtf b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(L1_gap L2_gap demo L1_stru $CTRL_FINAL) order(L1_gap L2_gap demo L1_stru $CTRL_FINAL) ///
    coeflabels(L1_gap "L1_gap" L2_gap "L2_gap" demo "demo" ///
      L1_stru "L1_stru" ln_pc_stu "ln_pc_stu" ///
      pc_retail_10k "pc_retail_10k" ln_mobile "ln_mobile" ///
      ln_pc_loan "ln_pc_loan") ///
    mtitles("正交两步System-GMM：统一1—2期") ///
    stats(N Groups Instruments DynSum DynJoint_p AR1_p AR2_p Hansen_p, ///
      labels("N" "Groups" "Instruments" "Dynamic_sum" "Dynamic_joint_p" "AR(1)_p" "AR(2)_p" "Hansen_p") fmt(0 0 0 4 4 4 4 4)) ///
    title("区域经济差异的部分动态System-GMM：统一1—2期工具变量") ///
    addnotes("样本为278个城市，估计期为2009—2024年；因变量与新基准回归一致。", ///
      "所有GMM-style变量（L1_gap、L2_gap、demo、L1_stru）统一使用1—2期滞后值并折叠。", ///
      "控制变量为ln_pc_stu、pc_retail_10k、ln_mobile和ln_pc_loan，与新基准回归一致。", ///
      "控制变量和年份虚拟变量作为外生工具变量；不保存DTA文件。")
log close
