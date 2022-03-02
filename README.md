**该流程将实现以下分析内容：**

##### PartI

收集人肝癌肿瘤组织(未治疗)单细胞测序分析数据(样本需包含肿瘤组织和免疫细胞)

1. QC所有样本，并剔除不合格数据
2. 分析并打分每个患者肿瘤样本中的CD3+T/CD8+T细胞的浸润比例(此处指T细胞占CD45+细胞比例)
3. 按照CD8+T(CD3+CD45+CD8+CD11b-B220-CD11C-NK1.1-gdT-)细胞浸润分数将不同患者分为两组：CD8+T high和CD8+T low组
   1. 针对CD8+T low组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   2. 针对CD8+T high组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   3. 在其他实体瘤scRNA数据库和TCGA data library中进行验证。
4. 按照CD3+T(CD3+CD45+CD11b-B220-CD11C-NK1.1-gdT-)细胞浸润分数将不同患者分为两组：CD3+T high和CD3+T low组
   1. 针对CD3+T low组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   2. 针对CD3+T high组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   3. 在其他实体瘤scRNA数据库和TCGA data library中进行验证。
5. 按照CD4+T(CD3+CD45+CD11b-B220-CD11C-NK1.1-gdT-)细胞浸润分数将不同患者分为两组：CD3+T high和CD3+T low组
   1. 针对CD4+T low组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   2. 针对CD4+T high组，寻找肿瘤组织中(CD45-, CAF, TEC, Malignant cell)特异性高表达/低表达的基因和信号通路；
   3. 在其他实体瘤scRNA数据库和TCGA data library中进行验证。

##### PartII

1. 针对CD3+T、CD4+T和CD8+T细胞中UP和Down基因(CAF、TEC、Malignant cell )，绘制基因表达热图/柱状图,（按照患者免疫浸润程度由高到低将Down和UP组的基因表达量进行分析）
2. 寻找CD3+T/CD4+T/CD8+T中CAF、TEC和Malignant cell中上调和下调的overlap进行，并进行表达热图展示和功能富集分析
3. CAF、TEC和Malignant cell分别在high组和low组中高表达基因的overlap分析，并进行功能富集分析

##### PartIII

1. 对差异基因进行TCGA的pan-cancer、survival和单细胞数据库验证
2. 对找到的差异基因进行与下列signatures的相关性分析：
   1. CD4: (1) CD3+, CD4+; (2) CD3+, CD4+, PD-1+, TIM3+, LAG3+; (3) CD3+, CD4+, IFNG+; (4) CD3+, CD4+, GranzymeB; (5) CD3+, CD4+, TNFG+;
   2. CD8: (1) CD3+, CD8+; (2) CD3+, CD8+, PD-1+, TIM3+, LAG3+; (3) CD3+, CD8+, IFNG+; (4) CD3+, CD8+, GranzymeB; (5) CD3+, CD8+, TNFA+
3. 找出CD8T(CD45+, CD3+, CD8+)中新细胞亚群，并进行差异基因功能分析
4. 找出CD4T(CD45+, CD3+, CD4+)中新细胞亚群，并进行差异基因功能分析
5. 将整个T细胞(CD45+, CD3+, CD4+, CD8+)划分为Tn、Tcm、Tem、Temra和Tex等

##### PartIV

1. 分析TAM中高表达的基因如LAIR1，并分析它们的表达在本数据集以及在TCGA肝癌数据库中与CD163，CD206, CD8，IFNg的相关性如下图（Human TAM：CD45+CD3e-CD19-CD56-CD16-CD206+CD163+ )
2. 分析G-MDSC中高表达的基因如X，并分析它们的表达在本数据集以及在TCGA数据库中与CD33，CD66b, CD11B, CD8， IFNg基因细胞相关性:（Human G-MDSC：CD45+ CD3-B220-NK1.1-gdTCR- CD11b+ CD14– CD15+ CD33+ CD66b+）
3. 分析M-MDSC中高表达的基因如X，并分析它们的表达在本数据集以及在TCGA数据库中与CD33，CD11B, CD8基因相关性:
   （Human M-MDSC：CD45+CD3-B220-NK1.1-gdTCR- CD11b+ CD14+ CD15– CD33+ CD66b–)
4. 分析Tex CD8+T中高表达的基因如X，并分析它们的表达在本数据集以及在TCGA数据库中与CD8，T-bet，IFNg,GranzymeB,IL-2, PD-1,LAG3,TIM3,TOX基因相关性 ( Tex CD8+T:CD45+CD8+CD3+CD152+ CD223+ CD244+ CD279+ CD366+ICOS+ TIGIT+ VISTA+Tox+)
5. 分析Tex CD4+T中高表达的基因如X，并分析它们的表达在本数据集以及在TCGA数据库中与CD4, T-bet，IFNg,GranzymeB,IL-2, PD-1,LAG3,TIM3,TOX基因相关性 ( Tex CD4+T:CD45+CD4+CD3+CD152+ CD223+ CD244+ CD279+ CD366+ICOS+ TIGIT+ VISTA+Tox+) ICOS+ TIGIT+ VISTA+Tox+)
6. 分析Treg中高表达的基因如X，并分析它们的表达在本数据集以及在TCGA数据库中与 FoxP3, Helios，CD8，CD11B ，CD33基因相关性( Treg:CD3+CD4+CD25+FoxP3+ cells )

##### PartV

1. 将所有HCC患者T细胞进行新的免疫分群10~16种亚型
   1. 在CD3+T高浸润和低浸润的患者中分析哪些新亚群富集在CD3+T高浸润肿瘤中，哪些新亚群富集在CD3+T低浸润肿瘤中；
   2. 对于在CD3+T高浸润/低浸润肿瘤中富集的新CD8+T细胞亚型，分析该亚型的特异性基因表达差异(相对其他15种T细胞），并分析该群细胞的特征基因与CTL，Teff, Tem，Tcm，Tex细胞的特征基因相似性;
   3. 对于在CD3+T高浸润/低浸润肿瘤中富集的新CD4+T细胞亚型，分析该亚型的特异性基因表达差异(相对其他15种T细胞），并分析该群细胞的特征基因与Th1，Th2，Treg，Th17，Tfh，和CD4+Teff, CD4+Tem，CD4+Tcm，CD4+Tex细胞的特征基因相似性;
2. 将所有HCC患者髓样细胞（CD45+CD11B+CD3-CD19-CD56-CD16-)进行新的免疫分群10~16种亚型
   1. 在CD3+T/CD8+T/CD4+T高浸润和低浸润的患者中分析哪些新亚群富集在CD3+T/CD8+T/CD4+T高浸润肿瘤中，哪些新亚群富集在CD3+T/CD8+T/CD4+T低浸润肿瘤中；
   2. 对于在在CD3+T/CD8+T/CD4+T高浸润/低浸润肿瘤中富集的新髓样细胞细胞亚型，分析该亚型的特异性基因表达差异(相对其他髓样细胞），并分析该群细胞的特征基因与G-MDSC,M-MDSC,TAM的相似性。

##### PartVI

1. **Objective:  针对HCC肿瘤组织、正常组织、健康人外周血中以下免疫细胞差异表达基因**
   1. Human TAM: CD45+CD3e-CD19-CD56-CD16-CD206+CD163+HLA-DR+CD11C+CD14+
   2. Human G-MDSC：CD45+ CD3-B220-NK1.1-gdTCR- HLA-DR– CD11b+ CD14– CD15+ CD33+ CD66b+
   3. Human M-MDSC：CD45+CD3-B220-NK1.1-gdTCR- HLA-DR–CD11b+ CD14+ CD15– CD33+ CD66b–)
   4. TEX exhausted CD8+T cells：CD45+CD8+CD3+CD223+ CD279+ CD366+ (TIM-3) Tox+
   5. TEX exhausted CD4+T cells：CD45+CD4+CD3+CD223+ CD279+ CD366+ (TIM-3) Tox+
   6. CTL: CD45+CD3+ CD8+ CD107a+ (LAMP-1) GranzymeB+ IFN-γ+
   7. Th1: CD45+CD3+ CD4+Tbet+ IFN-γ+CXCR3+
2. **Objective:  针对HCC肿瘤组织、正常组织、健康人外周血中以下免疫细胞相关性分析**
   1. 分析肿瘤患者相比于正常组织中的TAM中高表达的基因如LAIR1，并分析它们的表达在我们的数据库以及在TCGA 肝癌数据库中与CD163，CD206, CD8，IFNg的相关性如下图（Human TAM：CD45+CD3e-CD19-CD56-CD16-CD206+CD163+ )
   2. 分析肿瘤患者相比于正常组织中的G-MDSC中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与CD33，CD66b, CD11B, CD8， IFNg基因细胞相关性: （Human G-MDSC：CD45+ CD3-B220-NK1.1-gdTCR- CD11b+ CD14– CD15+ CD33+ CD66b+）
   3. 分析肿瘤患者相比于正常组织中的M-MDSC中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与CD33，CD11B, CD8基因相关性:（Human M-MDSC：CD45+CD3-B220-NK1.1-gdTCR- CD11b+ CD14+ CD15– CD33+ CD66b–)
   4. 分析肿瘤患者相比于正常组织中的Tex CD8+T中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与CD8，T-bet，IFNg,GranzymeB,IL-2, PD-1,LAG3,TIM3,TOX基因相关性 ( Tex CD8+T:CD45+CD8+CD3+CD152+ CD223+ CD244+ CD279+ CD366+ICOS+ TIGIT+ VISTA+Tox+)
   5. 分析肿瘤患者相比于正常组织中的Tex CD4+T中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与CD4, T-bet，IFNg,GranzymeB,IL-2, PD-1,LAG3,TIM3,TOX基因相关性 ( Tex CD4+T:CD45+CD4+CD3+CD152+ CD223+ CD244+ CD279+ CD366+ICOS+ TIGIT+ VISTA+Tox+) ICOS+ TIGIT+ VISTA+Tox+)
   6. 分析肿瘤患者相比于正常组织中的CTL中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与 CD45, CD3，CD8，CD107a ，GranzymeB，IFN-γ基因相关性
   7. 分析肿瘤患者相比于正常组织中的Th1中高表达的基因如X，并分析它们的表达在我们的数据库以及在TCGA数据库中与 CD45, CD3，CD4，Tbet ，CXCR3，IFN-γ基因相关性





