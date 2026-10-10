* ============================================================================
* 调节效应：最终保留两组短缺口径
* 方案1：人均零售短缺；方案2：移动电话不足。
* 口径、基准控制变量、固定效应和聚类标准误沿用原短缺口径补充方案。
* 仅输出两张RTF，不保存或改写DTA数据。
* ============================================================================
version 18.0
clear all
set more off
set maxvar 32767

global DATA   "D:\Desktop\区域经济差异实验（9.24之后）\机制检验\机制检验_共同样本_278城市.dta"
global OUTDIR "D:\Desktop\区域经济差异实验（9.24之后）\机制检验（最后）"
capture mkdir "$OUTDIR"

capture which reghdfe
if _rc {
    display as error "未找到 reghdfe；请先安装并检查 Stata ado 路径。"
    exit 199
}
capture which esttab
if _rc {
    display as error "未找到 esttab；请先安装 estout。"
    exit 199
}

use "$DATA", clear
sort id year
isid id year
xtset id year

* 与基准回归一致的五个控制变量。
global CTRL_FINAL ln_pc_stu ln_pc_lib ln_pc_hosp ln_retail ln_pc_mobile

* 样本期和主变量完整样本口径与原机制代码一致。
keep if inrange(year,2009,2024)
egen byte miss_main = rowmiss(abs_gap demo L1_stru $CTRL_FINAL)
keep if miss_main==0

* 人均零售规模：零售总额/估算人口×10000，再缩放到万元口径。
capture drop mobile_pc retail_pc retail_pc_10k z_mobile shortage ///
    demo_short stru_short
gen double mobile_pc = mobile/pop_proxy*10000 if mobile>=0 & pop_proxy>0
gen double retail_pc = retail/pop_proxy*10000 if retail>=0 & pop_proxy>0
gen double retail_pc_10k = retail_pc/10000

* 移动电话强度标准化；其不足指数越大表示移动电话基础越不足。
quietly summarize mobile_pc
gen double z_mobile = (mobile_pc-r(mean))/r(sd)

* 确保两组调节效应使用同一估计样本。
egen byte miss_mod = rowmiss(retail_pc_10k z_mobile)
keep if miss_mod==0
quietly count
display as result "共同估计样本观测值 N=" r(N)

capture program drop run_shortage
program define run_shortage
    syntax, MOD(varname) LABEL(string) FILE(string) TITLE(string)

    * 反向编码：数值越大，代表相对于样本平均水平越短缺/不足。
    quietly summarize `mod'
    gen double shortage = r(mean)-`mod'
    gen double demo_short = demo*shortage
    gen double stru_short = L1_stru*shortage

    quietly reghdfe abs_gap demo L1_stru shortage demo_short stru_short ///
        $CTRL_FINAL, absorb(id year) vce(cluster id)

    quietly test demo
    local pdemo = r(p)
    quietly test L1_stru
    local pstru = r(p)
    quietly test demo_short
    local pdemo_int = r(p)
    quietly test stru_short
    local pstru_int = r(p)
    estadd scalar Overall_R2=e(r2)
    estadd scalar P_demo=`pdemo'
    estadd scalar P_stru=`pstru'
    estadd scalar P_demo_interaction=`pdemo_int'
    estadd scalar P_stru_interaction=`pstru_int'

    esttab using "`file'", replace rtf b(4) se(4) ///
        star(* 0.10 ** 0.05 *** 0.01) ///
        keep(demo L1_stru shortage demo_short stru_short) ///
        order(demo L1_stru shortage demo_short stru_short) ///
        coeflabels(demo "demo" L1_stru "L1_stru" ///
            shortage "`label'短缺/不足指数" ///
            demo_short "demo×`label'短缺/不足" ///
            stru_short "L1_stru×`label'短缺/不足") ///
        mtitles("`label'短缺/不足") ///
        stats(N Overall_R2 P_demo P_stru P_demo_interaction P_stru_interaction, ///
            fmt(0 4 4 4 4 4) ///
            labels("观测值" "总体R²（reghdfe）" "demo p值" "L1_stru p值" ///
                "demo交互项p值" "L1_stru交互项p值")) ///
        title("`title'") ///
        addnotes("因变量abs_gap；控制五个基准控制变量、城市和年份固定效应，标准误按城市聚类。", ///
            "调节指数为样本均值减去原调节指标；值越大表示相对短缺/不足越严重。", ///
            "交互项刻画短缺/不足程度变化时demo和L1_stru边际效应的变化。", ///
            "规模/强度与短缺/不足为反向口径，不作为彼此独立的机制证据。")

    display as result "`label': demo交互项 p=" %7.4f `pdemo_int' ///
        "; L1_stru交互项 p=" %7.4f `pstru_int'
    drop shortage demo_short stru_short
end

* 方案1：人均零售短缺。
run_shortage, mod(retail_pc_10k) label("人均零售") ///
    file("$OUTDIR\方案1_人均零售短缺_连续交互.rtf") ///
    title("表A1·人均零售短缺调节效应")

* 方案2：移动电话不足。
run_shortage, mod(z_mobile) label("移动电话") ///
    file("$OUTDIR\方案2_移动电话不足_连续交互.rtf") ///
    title("表A2·移动电话不足调节效应")

display as result "两组短缺口径调节效应估计完成；仅输出两张RTF，未保存DTA。"
