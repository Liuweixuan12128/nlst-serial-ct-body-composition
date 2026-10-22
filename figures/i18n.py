"""Chinese labels for internal-review versions of the figures (FIG_LANG=zh). Submission figures stay English.

T(s) returns the Chinese string for an English template when FIG_LANG=zh, else s unchanged.
Templates may contain str.format placeholders; call T("... {n} ...").format(n=...).
"""
import os

LANG = os.environ.get("FIG_LANG", "en")

ZH = {
    # Figure 1
    "Public cohorts": "公共数据队列",
    "NLST, pT-staged   n = {n}": "NLST 病理分期   n = {n}",
    "   aged ≥65 y   n = {n}": "   其中≥65岁   n = {n}",
    "TCIA surgical   n = {n}": "TCIA 手术队列   n = {n}",
    "Scan–rescan pairs   n = {n}": "扫描-再扫描配对   n = {n}",
    "NSCLC resection": "非小细胞肺癌切除",
    "Outcomes": "结局",
    "Overall survival": "总生存",
    "Cancer vs non-cancer death": "肺癌 vs 非癌死亡",
    "Post-resection": "术后并发症",
    "complications": " ",
    "Lobar vs sublobar": "肺叶 vs 亚肺叶切除",
    "Routine chest CT": "常规胸部CT",
    "Deep-learning masks": "深度学习分割",
    "Vertebral-body levels": "椎体层面定位",
    "Tissue compartments": "组织分区",
    "Native (1.25 mm, lung kernel)": "原始图像（1.25 mm，肺算法）",
    "Harmonised": "噪声统一化后",
    "noise SD {sd} HU": "噪声SD {sd} HU",
    "fat\nwindow": "脂肪窗",
    "harmonised": "统一化后",
    "native": "原始",
    "Liver attenuation (HU)": "肝脏CT值（HU）",
    "Voxels (% per HU)": "体素（%/HU）",
    "Posterior 3D view": "后面观三维重建",
    "T5–T12 volume (red)": "T5–T12体积（红色）",
    "L": "左", "R": "右",
    "Muscle {a} cm², {b} HU": "肌肉 {a} cm²，{b} HU",
    "Skeletal muscle": "骨骼肌", "Intramuscular fat": "肌内脂肪", "Subcutaneous fat": "皮下脂肪",
    "Visceral/mediastinal fat": "内脏/纵隔脂肪",
    "T5–T12 volume\n{v} cm³": "T5–T12体积\n{v} cm³",
    "Distance above T12 centre (cm)": "距T12椎体中心高度（cm）",
    "Muscle area (cm²)": "肌肉面积（cm²）",
    "Patients from TCIA NSCLC-Radiogenomics (CC BY 3.0): a, c–e {a}; b {b}. Orange in b, voxels in the fat window "
    "(−190 to −30 HU).": "病例来自TCIA NSCLC-Radiogenomics（CC BY 3.0）：a、c–e为{a}；b为{b}。"
                        "b中橙色为落在脂肪窗（−190至−30 HU）内的体素。",
    # Figure 2
    "Staging CTs fully covering level (%)": "完整覆盖该层面的分期CT（%）",
    "n = {n} CTs": "n = {n} 例CT",
    "T5–T12 span": "T5–T12范围",
    "Soft kernel": "软组织算法", "Sharp kernel": "锐利算法",
    "{g}\nn = {n}": "{g}\nn = {n}",
    "Native": "原始", "Harm.": "统一化",
    "Image noise, liver SD (HU)": "图像噪声：肝脏SD（HU）",
    "Scan–rescan ICC (95% CI)": "扫描-再扫描ICC（95% CI）",
    "T5–T12 muscle volume": "T5–T12肌肉体积",
    "T12 muscle area": "T12肌肉面积",
    "T12 erector spinae area": "T12竖脊肌面积",
    "T5–T12 muscle density": "T5–T12肌肉密度",
    "T12 muscle density": "T12肌肉密度",
    "T12 intramuscular fat": "T12肌内脂肪",
    "T12 subcutaneous fat": "T12皮下脂肪",
    "T5–T12 subcutaneous fat": "T5–T12皮下脂肪",
    "T12 visceral fat": "T12内脏脂肪",
    "T5–T12 visceral/mediastinal fat": "T5–T12内脏/纵隔脂肪",
    "T12 trabecular bone density": "T12松质骨密度",
    "Mean of two CTs (cm³)": "两次CT均值（cm³）",
    "Scan 2 − scan 1 (cm³)": "第2次 − 第1次（cm³）",
    "bias {v}": "偏倚 {v}",
    "Scan 1 (HU)": "第1次CT（HU）", "Scan 2 (HU)": "第2次CT（HU）",
    "TCIA NSCLC-Radiogenomics. Scan 2: second staging CT of the same patient (mostly PET/CT), different protocol and "
    "kernel.": "TCIA NSCLC-Radiogenomics。第2次CT：同一患者的另一次分期CT（多为PET/CT），扫描协议和重建算法不同。",
    # Extended Data Fig. 1
    "Adjusted HR": "校正HR",
    "Overall death": "全因死亡",
    "Recurrence or death": "复发或死亡",
    "{t}\n({n} events)": "{t}\n（{n}个事件）",
    "Years after surgery": "术后年数",
    "Lowest quartile": "最低四分位", "Quartiles 2–4": "第2–4四分位",
    "T5–T12 muscle density": "T5–T12肌肉密度",
    "Muscle-to-vertebra ratio": "肌肉/椎体体积比",
    "T12 erector spinae density": "T12竖脊肌密度",
    "T12 intramuscular fat %": "T12肌内脂肪比例",
    "T12 bone density": "T12骨密度",
    "per SD lower": "每降低1 SD", "per SD higher": "每升高1 SD", ", ": "，",
    "TCIA NSCLC-Radiogenomics, age ≥65 y (n = {n}; {d} deaths). HRs adjusted for age, sex and pathological stage; "
    "filled squares, FDR q < 0.05. b, cut-off derived within this cohort (exploratory).":
        "TCIA NSCLC-Radiogenomics，≥65岁（n = {n}；死亡{d}例）。HR校正年龄、性别和病理分期；实心方块：FDR q < 0.05。"
        "b：切点在本队列内部确定（探索性）。",
    "Same contrast status": "增强状态相同", "Scan 1 contrast,\nscan 2 plain": "第1次增强，\n第2次平扫",
    "n = {n}; ICC {i}": "n = {n}；ICC {i}", "Density, scan 2 − scan 1 (HU)": "肌密度差：第2次 − 第1次（HU）",
    "Contrast explains part of the shift": "增强与否可解释部分偏差",
    "SYNTHETIC DATA — LAYOUT TEST ONLY": "合成数据 — 仅用于版式测试，不是结果",
    # public results figures (3-5)
    "{d} per decade": "每10岁 {d}", "Age at diagnosis (years)": "诊断时年龄（岁）", "Men": "男性", "Women": "女性",
    "Difference per decade of age (SD)": "每增加10岁的差异（SD）", "Stage": "分期", "z-score": "z分数",
    "NLST CT arm, pathologically staged lung cancer; last screening CT before diagnosis (n = {n}; {m} aged ≥65). "
    "Lines: local regression. b, adjusted for sex, kernel and current smoking; filled squares, FDR q < 0.05.":
        "NLST CT组、有病理分期的肺癌；诊断前最后一次筛查CT（n = {n}；≥65岁{m}例）。曲线为局部回归。"
        "b：校正性别、重建算法和当前吸烟；实心方块：FDR q < 0.05。",
    "lower": "越低", "higher": "越高", "{e}, SD\n(right = {w})": "{e}（SD）\n（向右{w}）",
    "Adjusted HR for death": "死亡校正HR", "All ages": "全年龄", "HR per SD": "每SD的HR",
    "By cohort, age ≥65 y": "分队列（≥65岁）",
    "TCIA NSCLC-Radiogenomics, TCGA-LUAD/LUSC and CPTAC (n = {n}; {m} aged ≥65). Cox models stratified by cohort, "
    "adjusted for age, sex, stage, histology, contrast and kernel; exposures z-scored within cohort and sex; filled "
    "squares, FDR q < 0.05.": "TCIA NSCLC-Radiogenomics、TCGA-LUAD/LUSC与CPTAC（n = {n}；≥65岁{m}例）。Cox模型按队列分层，"
                              "校正年龄、性别、分期、组织学、增强和重建算法；暴露按队列和性别标准化；实心方块：FDR q < 0.05。",
    "T12 erector spinae area": "T12竖脊肌面积", "Trunk cut by field of view (%)": "躯干被视野截断（%）",
    "Muscle-to-vertebra ratio*": "肌肉/椎体体积比*", "T12 subcutaneous fat*": "T12皮下脂肪*",
    "ESM area": "竖脊肌面积", "ESM density": "竖脊肌密度",
    "TCIA NSCLC-Radiogenomics, TCGA-LUAD/LUSC and CPTAC (n = {n}; {m} aged ≥65). Cox models stratified by cohort, "
    "adjusted for age, sex, stage, histology, contrast and kernel; exposures z-scored within cohort and sex; filled "
    "squares, FDR q < 0.05. *CTs without field-of-view truncation.": "TCIA NSCLC-Radiogenomics、TCGA-LUAD/LUSC与CPTAC"
    "（n = {n}；≥65岁{m}例）。Cox模型按队列分层，校正年龄、性别、分期、组织学、增强和重建算法；暴露按队列和性别标准化；"
    "实心方块：FDR q < 0.05。*仅视野未截断的CT。",
    "Each cohort held out in turn; models fitted on the others. Clinical: age, sex, stage, histology. Body composition: "
    "erector spinae area and density, intramuscular fat %.": "每个队列轮流留出作测试，模型在其余队列中拟合。临床变量："
    "年龄、性别、分期、组织学。体成分：竖脊肌面积和密度、肌内脂肪比例。",
    "Age ≥65 y": "≥65岁",
    "deaths": "例死亡", "Muscle ratio": "肌肉比", "Muscle density": "肌密度", "IMAT %": "肌内脂肪%", "Harrell's C in the held-out cohort": "留出队列中的Harrell C指数", "Clinical model": "临床模型",
    "Clinical + body composition": "临床+体成分", "ΔC (95% bootstrap CI)": "ΔC（95% bootstrap CI）",
    "3-year time-dependent AUC": "3年时间依赖AUC",
    "Each cohort held out in turn; models fitted on the others. Clinical: age, sex, stage, histology. Body composition: "
    "muscle-to-vertebra ratio, muscle density, intramuscular fat %.": "每个队列轮流留出作测试，模型在其余队列中拟合。"
                                                                     "临床变量：年龄、性别、分期、组织学。体成分：肌肉/椎体体积比、肌肉密度、肌内脂肪比例。",
    # Figures 3-6 (templates)
    "3D muscle index": "三维肌指数", "T12 muscle index": "T12肌指数", "Pectoralis index (T4)": "胸大肌指数（T4）",
    "Erector spinae index (T12)": "竖脊肌指数（T12）", "Intramuscular fat (T12)": "肌内脂肪（T12）",
    "Visceral fat (T12)": "内脏脂肪（T12）", "Subcutaneous fat (T12)": "皮下脂肪（T12）",
    "Vertebral bone density (T12)": "椎体骨密度（T12）",
    "Clavien–Dindo ≥III": "Clavien–Dindo ≥III级", "Pulmonary complications": "肺部并发症", "90-day death": "90天死亡",
    "Prolonged stay": "住院时间延长",
    "Normal": "正常", "Low muscle only": "仅低肌量", "Low density only": "仅低肌密度",
    "Low muscle and density": "低肌量+低肌密度",
    "All patients": "全部患者", "Low muscle": "低肌量", "Normal muscle": "正常肌量", "Low muscle density": "低肌密度",
    "Normal muscle density": "正常肌密度", "Age ≥75 y": "年龄≥75岁", "Age 65–74 y": "年龄65–74岁",
    "IA ≤2 cm, low muscle": "IA期≤2 cm，低肌量", "IA ≤2 cm, normal muscle": "IA期≤2 cm，正常肌量",
    "Development": "开发队列", "Temporal validation": "时间验证", "External validation": "外部验证",
    "{e}\n(sex-specific z-score)": "{e}\n（按性别标准化z分数）",
    "Adjusted OR, Clavien–Dindo ≥III": "校正OR（Clavien–Dindo ≥III级）",
    "Adjusted OR": "校正OR", "Overall survival": "总生存", "Cumulative incidence": "累积发生率",
    "Non-cancer death": "非癌死亡", "Cancer death": "癌死亡", "Number at risk": "风险人数",
    "Subdistribution HR": "亚分布HR", "Cause-specific death": "按死因分类的死亡",
    "Non-cancer": "非癌", "Cancer": "癌",
    "3D muscle index, per SD lower": "三维肌指数，每降低1 SD", "T12 muscle density, per SD lower": "T12肌密度，每降低1 SD",
    "1 − specificity": "1 − 特异度", "Sensitivity": "敏感度", "Predicted risk": "预测风险", "Observed risk": "实际风险",
    "Threshold probability": "阈值概率", "Net benefit": "净获益",
    "Clinical model (M0)": "临床模型（M0）", "Clinical + body composition (M1)": "临床+体成分（M1）",
    "Treat all": "全部干预", "Treat none": "全不干预",
    "AUC / Harrell's C": "AUC / Harrell C指数",
    "Overall survival (IPTW)": "总生存（IPTW加权）", "Lobar": "肺叶切除", "Sublobar": "亚肺叶切除",
    "HR sublobar vs lobar": "HR（亚肺叶 vs 肺叶）", "IPTW HR {h} ({l}–{u})": "IPTW HR {h}（{l}–{u}）",
    "{t} (n = {n})": "{t}（n = {n}）",
    "{c}\nn = {n}; {e} events": "{c}\nn = {n}；{e}个事件",
    "AUC {a} vs {b}": "AUC {a} vs {b}",
    "CD ≥III, {c}": "CD ≥III，{c}", "Survival, {c}": "生存，{c}",
    "development": "开发", "temporal": "时间验证", "external": "外部验证",
    "Death": "死亡",
    "P nonlinear": "非线性P", "P interaction, low muscle": "交互作用P（低肌量）",
    "Filled squares, FDR q < 0.05. Models adjusted for age, sex, BMI, Charlson index, ASA ≥III, ppoFEV1, approach, "
    "extent of resection and centre.": "实心方块：FDR q < 0.05。模型校正年龄、性别、BMI、Charlson指数、ASA ≥III、ppoFEV1、手术入路、切除范围和中心。",
    "Low muscle / low density: sex-specific lowest quartile in the development cohort. Filled squares, FDR q < 0.05.\n"
    "HRs adjusted for age, sex, BMI, Charlson, ASA, ppoFEV1, approach, resection, stage and adjuvant therapy, "
    "stratified by centre.": "低肌量/低肌密度：开发队列内按性别的最低四分位。实心方块：FDR q < 0.05。\n"
                             "HR校正年龄、性别、BMI、Charlson、ASA、ppoFEV1、入路、切除范围、分期和辅助治疗，并按中心分层。",
    "Stabilised inverse-probability-of-treatment weights; filled squares, 95% CI excludes 1.":
        "稳定化逆概率处理加权（IPTW）；实心方块：95% CI不含1。",
}



def T(s):
    return ZH.get(s, s) if LANG == "zh" else s
