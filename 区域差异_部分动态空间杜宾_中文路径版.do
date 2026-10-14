* 研究对象：278个城市，2009—2024年；因变量越小表示增长率差距越小。
* 模型：abs_gap(t)=rho*W*abs_gap(t)+tau*abs_gap(t-1)
*       +beta_d*demo(t)+beta_s*stru(t-1)+theta_s*W*stru(t-1)
*       +四项控制变量+城市固定效应+省份×年份固定效应+误差。
* 仅给结构变量设置空间杜宾项，demo及控制变量没有W项。
* 文献依据：Barabuffi等(2025)反距离平方权重；Zhang和Gibson(2025)
* 空间计量与直接/间接效应；Ascani等(2025)只对关键变量设置空间滞后。
* 引用详见文末。xsmle结果与前面的Python SPJ-QML不是同一估计器。
* 这些文献支持方法选择，不保证本研究的效应方向或显著性。
version 18.0
clear all
set more off

* 若原始数据文件夹移动，只修改ROOT；OUT为完整和精简结果的输出目录。
global ROOT "D:/Desktop/工作/研究生师门工作/专利转移/【最终成功】专利转移与区域经济增长【2009-2024年】"
global OUT "$ROOT/空间SDM_仅结构项"
capture mkdir "$OUT"
capture log close _all
log using "$OUT/运行结果.log", text replace

capture which xsmle
if _rc {
    display as error "缺少 xsmle。请先在 Stata 中运行：ssc install xsmle"
    exit 199
}
capture which esttab
if _rc {
    display as error "缺少 esttab。请先在 Stata 中运行：ssc install estout"
    exit 199
}

* 1. 因变量：城市实际人均GDP增长率与当年城市均值的绝对差。
* 必须先在原始2005—2024数据中算年度均值，再筛278城；顺序与基准一致。
* L1_stru是上一年结构增长；控制变量沿用最终基准回归的四项。
use "$ROOT/删除部分城市和缩尾1%的数据.dta", clear
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

merge m:1 id using "$ROOT/省份城市id对照.dta", ///
    keepusing(province) keep(master match) nogen
encode province, gen(province_id)
local controls ln_pc_stu pc_retail_10k ln_mobile ln_pc_loan
* 只保留2009—2024年每年都满足基准变量完整性的城市。
* 2008年真实观测留作L.abs_gap的支持期，不作为目标估计年份。
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

* 2. 原始单位系数很小：demo和L1_stru分别除以基准样本标准差。
* 新系数对应解释变量增加1个标准差；标准误同比放大，显著性不变。
* 标准差只用2009—2024年4448条，不用2008支持期。
quietly summarize demo if inrange(year,2009,2024)
scalar sd_demo=r(sd)
quietly summarize L1_stru if inrange(year,2009,2024)
scalar sd_stru=r(sd)
assert scalar(sd_demo)>0 & scalar(sd_stru)>0
gen double demo_sd=demo/scalar(sd_demo)
gen double stru_sd=L1_stru/scalar(sd_stru)
display as result "demo 标准差 = " %12.6f scalar(sd_demo)
display as result "L1_stru 标准差 = " %12.6f scalar(sd_stru)

* 3. 显式生成省份×年份虚拟变量，吸收各省逐年的共同冲击。
* 城市固定效应由fe type(ind)吸收；2009年是省内基准年。
quietly summarize province_id, meanonly
local np=r(max)
local py ""
forvalues p=1/`np' {
    forvalues t=2010/2024 {
        gen byte py_`p'_`t'=(province_id==`p' & year==`t')
        local py "`py' py_`p'_`t'"
    }
}

* 4. 空间权重：用城市经纬度求球面公里距离，非自身W_ij=1/d_ij^2。
* 对角线为0，各行除以行和，故每行和为1；城市顺序按id。
* 反距离平方参考Barabuffi等(2025)，此处城市单位和278城样本是本研究设定。
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

* 5. dlag(1)加入上一期abs_gap；wmat()加入同期W*abs_gap。
* durbin(stru_sd)只加入W*L1_stru，不加入W*demo或W*控制变量。
* effects输出短期和长期效应；只看回归系数不足以判断总作用。
set seed 20261014
xsmle abs_gap demo_sd stru_sd `controls' `py', ///
    wmat(W_INV2) model(sdm) fe type(ind) dlag(1) ///
    durbin(stru_sd) effects vceeffects(sim, nsim(500)) ///
    nolog iterate(500)

count if e(sample)
if r(N)!=4448 {
    display as error "有效样本不是 4448；请查看年份及 xsmle 版本，不要直接引用结果。"
    tab year if e(sample)
    exit 459
}
tab year if e(sample)
matrix list e(b)

* 效应解读：短期=[I-rho*W]^-1*(beta*I+theta*W)。
* 长期=[(1-tau)*I-rho*W]^-1*(beta*I+theta*W)；demo的theta=0。
* 直接效应为平均对角元，总效应为平均行和，间接效应=总效应-直接效应。
* 负的总效应代表abs_gap降低；统计显著性仍须看标准误和p值。

* 6. 完整输出：原生估计文件、所有参数与效应的 RTF/CSV，以及运行日志。
estimates store SDM_FULL
estimates save "$OUT/完整估计.ster", replace
esttab SDM_FULL using "$OUT/完整结果.rtf", replace rtf ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N, fmt(0) labels("估计样本量")) ///
    title("部分动态空间杜宾模型：全部系数与效应")
esttab SDM_FULL using "$OUT/完整结果.csv", replace csv ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01)

* 7. 精简表只列四项：两变量各自的短期总效应和长期总效应。
* 完整表还包含所有回归系数、直接/间接效应、标准误和p值。
local core SR_Total:demo_sd SR_Total:stru_sd ///
    LR_Total:demo_sd LR_Total:stru_sd
esttab SDM_FULL using "$OUT/核心短长期总效应.rtf", replace rtf ///
    keep(`core') order(`core') ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    coeflabels(SR_Total:demo_sd "demo：短期总效应" ///
               SR_Total:stru_sd "L1_stru：短期总效应" ///
               LR_Total:demo_sd "demo：长期总效应" ///
               LR_Total:stru_sd "L1_stru：长期总效应") ///
    stats(N, fmt(0) labels("估计样本量")) ///
    title("核心变量的短期与长期总效应（一标准差变化）")
esttab SDM_FULL using "$OUT/核心短长期总效应.csv", replace csv ///
    keep(`core') order(`core') ///
    cells("b(star fmt(6)) se(par fmt(6)) p(par fmt(4))") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    coeflabels(SR_Total:demo_sd "demo：短期总效应" ///
               SR_Total:stru_sd "L1_stru：短期总效应" ///
               LR_Total:demo_sd "demo：长期总效应" ///
               LR_Total:stru_sd "L1_stru：长期总效应")
display as result "输出目录：$OUT"
log close

* 主要参考文献（DOI可用于查原文；页码指用户上传PDF的页序）：
* Barabuffi S, Cricchio J, Di Minin A (2025). The 'picking the fittest'
* approach and spatial dynamics in China's artificial intelligence regional
* development. Papers in Regional Science 104, 100096.
* DOI: 10.1016/j.pirs.2025.100096；反距离平方矩阵见上传2.pdf第6页。
* Zhang X, Gibson J (2025). Local economic effects of connecting to China's
* high-speed rail network: evidence from spatial econometric models.
* The Annals of Regional Science 74, article 56.
* DOI: 10.1007/s00168-025-01384-4；空间效应分解见上传6.pdf第14—15页。
* Ascani A, Bettarelli L, Resmini L (2025). Unpacking the black box of
* knowledge externalities: intra-MNE networks vs. spatial spillovers.
* Papers in Regional Science 104, 100115.
* DOI: 10.1016/j.pirs.2025.100115；重点解释变量空间项见上传7.pdf第6—8页。
* LeSage JP, Pace RK (2009). Introduction to Spatial Econometrics.
* Chapman and Hall/CRC；直接、间接、总效应的经典定义。
* Belotti F, Hughes G, Piano Mortari A (2017). Spatial panel-data models
* using Stata. The Stata Journal 17(1):139—180；xsmle命令的方法来源。
* 本脚本另沿用用户上传的空间.do中xsmle程序结构，但样本和变量按本项目重建。
