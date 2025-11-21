# Lab 14 – Markdown Version

BIMM 143 – Lab 15: Pathway Analysis from RNA-Seq Results
Lulu Zha (PID: a17942879)
2025-11-19
Contents
Submission Header 2
Overview 2
Section 1 – Differential Expression Analysis (DESeq2) 2
Load libraries and data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 2
Remove the length column . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 3
Filter out genes with all zeros . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 3
Create DESeqDataSet and run DESeq2 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4
Get results and summarize . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4
Section 1 – Volcano Plot 5
Section 1 – Adding Gene Annotation 7
Save ordered results to CSV . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 9
Section 2 – KEGG Pathway Analysis with GAGE and pathview 9
Prepare foldchange vector for GAGE . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 10
Run GAGE . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 10
Pathview example for a down-regulated pathway . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
Top 5 up-regulated pathways (greater) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
Q7 – Top 5 down-regulated pathways with pathview . . . . . . . . . . . . . . . . . . . . . . . . . . 11
Section 3 – Gene Ontology (GO) with GAGE 12
Section 4 – Reactome Analysis 13
Section 5 – GO Online Enrichment (Optional) 14
Session Information (Reproducibility) 15
1
Submission Header
Name:Lulu Zha
PID:a17942879
Course:BIMM 143 – Bioinformatics Lab
Lab:Lecture 15 – Pathway Analysis from RNA-Seq Results
Overview
In this lab I practice taking RNA-seq differential expression results and using them to dopathway analysis.
Big picture steps:
•Import raw count data and sample metadata.
•Run DESeq2 to find differentially expressed genes.
•Make a volcano plot and annotate genes (Ensembl→SYMBOL / Entrez / gene name).
•RunGAGEfor KEGG pathway enrichment and plot pathways withpathview.
•Do additional pathway / GO analysis usingReactomeand the online GO enrichment tools.
This report includes the R code I used plus short answers to each question in the lab handout.
Section 1 – Differential Expression Analysis (DESeq2)
In this section I load the count matrix and metadata, clean the data (remove the length column and all-zero
rows), run DESeq2, and look at the summary.
Load libraries and data
library(DESeq2)
I keep the CSV files in thesame folderas this Rmd file, with the names that show up in Finder
(GSE37704_metadata (1).csvandGSE37704_featurecounts (1).csv).
metaFile <- "GSE37704_metadata (1).csv"
countFile <- "GSE37704_featurecounts (1).csv"
# Import metadata and peek
colData <-read.csv(metaFile, row.names = 1)
head(colData)
## condition
## SRR493366 control_sirna
## SRR493367 control_sirna
2
## SRR493368 control_sirna
## SRR493369 hoxa1_kd
## SRR493370 hoxa1_kd
## SRR493371 hoxa1_kd
# Import count data
countData <-read.csv(countFile, row.names = 1)
head(countData)
## length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000186092 918 0 0 0 0 0
## ENSG00000279928 718 0 0 0 0 0
## ENSG00000279457 1982 23 28 29 29 28
## ENSG00000278566 939 0 0 0 0 0
## ENSG00000273547 939 0 0 0 0 0
## ENSG00000187634 3214 124 123 205 207 212
## SRR493371
## ENSG00000186092 0
## ENSG00000279928 0
## ENSG00000279457 46
## ENSG00000278566 0
## ENSG00000273547 0
## ENSG00000187634 258
Remove the length column
Q1 – Complete the code to remove the first (length) column fromcountData:
# Remove the odd first $length column
countData <-as.matrix(countData[,-1])
head(countData)
## SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
## ENSG00000186092 0 0 0 0 0 0
## ENSG00000279928 0 0 0 0 0 0
## ENSG00000279457 23 28 29 29 28 46
## ENSG00000278566 0 0 0 0 0 0
## ENSG00000273547 0 0 0 0 0 0
## ENSG00000187634 124 123 205 207 212 258
Answer (Q1):I usecountData[, -1]to drop the first column (the length column) and keep only the
count columns.
Filter out genes with all zeros
Q2 – FiltercountDatato exclude genes (rows) where all samples have 0 reads:
# Filter out rows (genes) where the sum of counts across samples is 0
countData <- countData[rowSums(countData)>0, ]
head(countData)
3
## SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
## ENSG00000279457 23 28 29 29 28 46
## ENSG00000187634 124 123 205 207 212 258
## ENSG00000188976 1637 1831 2383 1226 1326 1504
## ENSG00000187961 120 153 180 236 255 357
## ENSG00000187583 24 48 65 44 48 64
## ENSG00000187642 4 9 16 14 16 16
Answer (Q2):rowSums(countData) > 0returns TRUE only for genes that have at least one non-zero
count across all samples. Using that as a row filter keeps only genes with some data.
Create DESeqDataSet and run DESeq2
dds <-DESeqDataSetFromMatrix(countData = countData,
colData = colData,
design =~condition)
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors
dds <-DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
dds
## class: DESeqDataSet
## dim: 15975 6
## metadata(1): version
## assays(4): counts mu H cooks
## rownames(15975): ENSG00000279457 ENSG00000187634 ... ENSG00000276345
## ENSG00000271254
## rowData names(22): baseMean baseVar ... deviance maxCooks
## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
## colData names(2): condition sizeFactor
Get results and summarize
4
res <-results(dds, contrast =c("condition", "hoxa1_kd", "control_sirna"))
Q3 – Callsummary()on the results:
summary(res)
##
## out of 15975 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up) : 4349, 27%
## LFC < 0 (down) : 4396, 28%
## outliers [1] : 0, 0%
## low counts [2] : 1237, 7.7%
## (mean count < 0)
## [1] see ’cooksCutoff’ argument of ?results
## [2] see ’independentFiltering’ argument of ?results
Answer (Q3):I runsummary(res). This tells me how many genes are significantly up-regulated and
down-regulated at the default adjusted p-value cutoff (0.1), plus how many are low count or outliers.
Section 1 – Volcano Plot
Here I make a volcano plot of log2 fold change vs –log(p-adj) and then improve the aesthetics with color.
plot(res$log2FoldChange,-log(res$padj))
5
−4 −2 0 2 4 6 80200400600
res$log2FoldChange−log(res$padj)Q4 – Improve the plot with colors and axis labels:
# Color vector for all genes
mycols <-rep("gray",nrow(res))
# Color red for genes with |log2FC| > 2
mycols[abs(res$log2FoldChange)>2] <- "red"
# Color blue for genes with padj < 0.01 and |log2FC| > 2
inds <- (res$padj<0.01)&(abs(res$log2FoldChange)>2)
mycols[inds] <- "blue"
plot(res$log2FoldChange,
-log(res$padj),
col = mycols,
xlab = "Log2(FoldChange)",
ylab = "-Log(P-value)")
6
−4 −2 0 2 4 6 80200400600
Log2(FoldChange)−Log(P−value)Answer (Q4):- The logical condition insideindsisres$padj < 0.01 & abs(res$log2FoldChange) >
2. - I usecol = mycolsinsideplot()so that strong hits are highlighted.
Section 1 – Adding Gene Annotation
Here I addSYMBOL,ENTREZID, andGENENAMEto my DESeq2 results usingmapIds()and the
org.Hs.eg.dbannotation.
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
## [1] "ACCNUM" "ALIAS" "ENSEMBL" "ENSEMBLPROT" "ENSEMBLTRANS"
## [6] "ENTREZID" "ENZYME" "EVIDENCE" "EVIDENCEALL" "GENENAME"
## [11] "GENETYPE" "GO" "GOALL" "IPI" "MAP"
## [16] "OMIM" "ONTOLOGY" "ONTOLOGYALL" "PATH" "PFAM"
## [21] "PMID" "PROSITE" "REFSEQ" "SYMBOL" "UCSCKG"
## [26] "UNIPROT"
Q5 – Fill in themapIds()calls:
7
# Add gene symbol
res$symbol <-mapIds(org.Hs.eg.db,
keys =row.names(res),
keytype = "ENSEMBL",
column = "SYMBOL",
multiVals = "first")
## ’select()’ returned 1:many mapping between keys and columns
# Add Entrez ID
res$entrez <-mapIds(org.Hs.eg.db,
keys =row.names(res),
keytype = "ENSEMBL",
column = "ENTREZID",
multiVals = "first")
## ’select()’ returned 1:many mapping between keys and columns
# Add full gene name
res$name <-mapIds(org.Hs.eg.db,
keys =row.names(res),
keytype = "ENSEMBL",
column = "GENENAME",
multiVals = "first")
## ’select()’ returned 1:many mapping between keys and columns
head(res, 10)
## log2 fold change (MLE): condition hoxa1_kd vs control_sirna
## Wald test p-value: condition hoxa1 kd vs control sirna
## DataFrame with 10 rows and 9 columns
## baseMean log2FoldChange lfcSE stat pvalue
## <numeric> <numeric> <numeric> <numeric> <numeric>
## ENSG00000279457 29.913579 0.1792571 0.3248215 0.551863 5.81042e-01
## ENSG00000187634 183.229650 0.4264571 0.1402658 3.040350 2.36304e-03
## ENSG00000188976 1651.188076 -0.6927205 0.0548465 -12.630156 1.43993e-36
## ENSG00000187961 209.637938 0.7297556 0.1318599 5.534326 3.12428e-08
## ENSG00000187583 47.255123 0.0405765 0.2718928 0.149237 8.81366e-01
## ENSG00000187642 11.979750 0.5428105 0.5215598 1.040744 2.97994e-01
## ENSG00000188290 108.922128 2.0570638 0.1969053 10.446970 1.51281e-25
## ENSG00000187608 350.716868 0.2573837 0.1027266 2.505522 1.22271e-02
## ENSG00000188157 9128.439422 0.3899088 0.0467164 8.346302 7.04333e-17
## ENSG00000237330 0.158192 0.7859552 4.0804729 0.192614 8.47261e-01
## padj symbol entrez name
## <numeric> <character> <character> <character>
## ENSG00000279457 6.86555e-01 NA NA NA
## ENSG00000187634 5.15718e-03 SAMD11 148398 sterile alpha motif ..
## ENSG00000188976 1.76553e-35 NOC2L 26155 NOC2 like nucleolar ..
## ENSG00000187961 1.13413e-07 KLHL17 339451 kelch like family me..
## ENSG00000187583 9.19031e-01 PLEKHN1 84069 pleckstrin homology ..
8
## ENSG00000187642 4.03379e-01 PERM1 84808 PPARGC1 and ESRR ind..
## ENSG00000188290 1.30538e-24 HES4 57801 hes family bHLH tran..
## ENSG00000187608 2.37452e-02 ISG15 9636 ISG15 ubiquitin like..
## ENSG00000188157 4.21970e-16 AGRN 375790 agrin
## ENSG00000237330 NA RNF223 401934 ring finger protein ..
Answer (Q5):- Forsymbol:keys = row.names(res),keytype = "ENSEMBL",column = "SYMBOL". -
Forentrez:keys = row.names(res),column = "ENTREZID". - Forname:keytype = "ENSEMBL",column
= "GENENAME".
Save ordered results to CSV
Q6 – Reorder by p-value and write to CSV:
# Order the results by p-value (as shown in the lab text)
res <- res[order(res$pvalue), ]
# Save as a CSV file
write.csv(as.data.frame(res), file = "deseq_results.csv")
Answer (Q6):I reorderresbyres$pvalueand then use
write.csv(as.data.frame(res), file = "deseq_results.csv")
to save the full annotated table in my project directory.
Section 2 – KEGG Pathway Analysis with GAGE and pathview
In this section I useGAGEto find KEGG pathways enriched for up- or down-regulated genes and then
pathviewto plot those pathways with my RNA-seq fold changes overlaid.
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
# Focus on signaling and metabolic pathways only
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
## $‘hsa00232 Caffeine metabolism‘
## [1] "10" "1544" "1548" "1549" "1553" "7498" "9"
##
## $‘hsa00983 Drug metabolism - other enzymes‘
## [1] "10" "1066" "10720" "10941" "151531" "1548" "1549" "1551"
## [9] "1553" "1576" "1577" "1806" "1807" "1890" "221223" "2990"
## [17] "3251" "3614" "3615" "3704" "51733" "54490" "54575" "54576"
9
## [25] "54577" "54578" "54579" "54600" "54657" "54658" "54659" "54963"
## [33] "574537" "64816" "7083" "7084" "7172" "7363" "7364" "7365"
## [41] "7366" "7367" "7371" "7372" "7378" "7498" "79799" "83549"
## [49] "8824" "8833" "9" "978"
##
## $‘hsa00230 Purine metabolism‘
## [1] "100" "10201" "10606" "10621" "10622" "10623" "107" "10714"
## [9] "108" "10846" "109" "111" "11128" "11164" "112" "113"
## [17] "114" "115" "122481" "122622" "124583" "132" "158" "159"
## [25] "1633" "171568" "1716" "196883" "203" "204" "205" "221823"
## [33] "2272" "22978" "23649" "246721" "25885" "2618" "26289" "270"
## [41] "271" "27115" "272" "2766" "2977" "2982" "2983" "2984"
## [49] "2986" "2987" "29922" "3000" "30833" "30834" "318" "3251"
## [57] "353" "3614" "3615" "3704" "377841" "471" "4830" "4831"
## [65] "4832" "4833" "4860" "4881" "4882" "4907" "50484" "50940"
## [73] "51082" "51251" "51292" "5136" "5137" "5138" "5139" "5140"
## [81] "5141" "5142" "5143" "5144" "5145" "5146" "5147" "5148"
## [89] "5149" "5150" "5151" "5152" "5153" "5158" "5167" "5169"
## [97] "51728" "5198" "5236" "5313" "5315" "53343" "54107" "5422"
## [105] "5424" "5425" "5426" "5427" "5430" "5431" "5432" "5433"
## [113] "5434" "5435" "5436" "5437" "5438" "5439" "5440" "5441"
## [121] "5471" "548644" "55276" "5557" "5558" "55703" "55811" "55821"
## [129] "5631" "5634" "56655" "56953" "56985" "57804" "58497" "6240"
## [137] "6241" "64425" "646625" "654364" "661" "7498" "8382" "84172"
## [145] "84265" "84284" "84618" "8622" "8654" "87178" "8833" "9060"
## [153] "9061" "93034" "953" "9533" "954" "955" "956" "957"
## [161] "9583" "9615"
Prepare foldchange vector for GAGE
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
head(foldchanges)
## 1266 54855 1465 2034 2150 6659
## -2.422719 3.201955 -2.313738 -1.888019 3.344508 2.392288
Run GAGE
keggres <-gage(foldchanges, gsets = kegg.sets.hs)
attributes(keggres)
## $names
## [1] "greater" "less" "stats"
# Look at first few down (less) pathways
head(keggres$less)
10
## p.geomean stat.mean p.val
## hsa04110 Cell cycle 8.995727e-06 -4.378644 8.995727e-06
## hsa03030 DNA replication 9.424076e-05 -3.951803 9.424076e-05
## hsa03013 RNA transport 1.375901e-03 -3.028500 1.375901e-03
## hsa03440 Homologous recombination 3.066756e-03 -2.852899 3.066756e-03
## hsa04114 Oocyte meiosis 3.784520e-03 -2.698128 3.784520e-03
## hsa00010 Glycolysis / Gluconeogenesis 8.961413e-03 -2.405398 8.961413e-03
## q.val set.size exp1
## hsa04110 Cell cycle 0.001448312 121 8.995727e-06
## hsa03030 DNA replication 0.007586381 36 9.424076e-05
## hsa03013 RNA transport 0.073840037 144 1.375901e-03
## hsa03440 Homologous recombination 0.121861535 28 3.066756e-03
## hsa04114 Oocyte meiosis 0.121861535 102 3.784520e-03
## hsa00010 Glycolysis / Gluconeogenesis 0.212222694 53 8.961413e-03
The top “less” (down-regulated) pathway isCell cycle (hsa04110)in this dataset.
Pathview example for a down-regulated pathway
# Example: Cell cycle pathway
pathview(gene.data = foldchanges,
pathway.id = "hsa04110",
species = "hsa")
This draws the KEGG cell cycle pathway and colors each gene node by its log2 fold-change.
Top 5 up-regulated pathways (greater)
# Top 5 upregulated pathways
keggrespathways <-rownames(keggres$greater)[1:5]
# Extract the 8-character KEGG IDs
keggresids <-substr(keggrespathways, start = 1, stop = 8)
keggresids
# Plot all 5 pathways (up) – can take some time
pathview(gene.data = foldchanges,
pathway.id = keggresids,
species = "hsa")
Q7 – Top 5 down-regulated pathways with pathview
Q7 – Do the same procedure to plot the pathview figures for the top 5 down-regulated path-
ways:
# Get the top 5 downregulated (less) pathways
keggrespathways.down <-rownames(keggres$less)[1:5]
11
# Extract the KEGG IDs (first 8 characters)
keggresids.down <-substr(keggrespathways.down, start = 1, stop = 8)
keggresids.down
# Draw pathview plots for these downregulated pathways
pathview(gene.data = foldchanges,
pathway.id = keggresids.down,
species = "hsa")
Answer (Q7):I repeat the same steps as for the “greater” set but usingkeggres$less.
- First get the top 5 pathways withrownames(keggres$less)[1:5].
- Pull out the KEGG IDs usingsubstr(..., 1, 8).
- Callpathview()withpathway.id = keggresids.downandspecies = "hsa".
Section 3 – Gene Ontology (GO) with GAGE
Now I run a similar enrichment analysis using GOBiological Process (BP)gene sets.
data(go.sets.hs)
data(go.subs.hs)
# Biological Process GO sets
gobpsets <- go.sets.hs[go.subs.hs$BP]
gobpres <-gage(foldchanges, gsets = gobpsets, same.dir = TRUE)
lapply(gobpres, head)
## $greater
## p.geomean stat.mean p.val
## GO:0007156 homophilic cell adhesion 8.519724e-05 3.824205 8.519724e-05
## GO:0002009 morphogenesis of an epithelium 1.396681e-04 3.653886 1.396681e-04
## GO:0048729 tissue morphogenesis 1.432451e-04 3.643242 1.432451e-04
## GO:0007610 behavior 1.925222e-04 3.565432 1.925222e-04
## GO:0060562 epithelial tube morphogenesis 5.932837e-04 3.261376 5.932837e-04
## GO:0035295 tube development 5.953254e-04 3.253665 5.953254e-04
## q.val set.size exp1
## GO:0007156 homophilic cell adhesion 0.1951953 113 8.519724e-05
## GO:0002009 morphogenesis of an epithelium 0.1951953 339 1.396681e-04
## GO:0048729 tissue morphogenesis 0.1951953 424 1.432451e-04
## GO:0007610 behavior 0.1967577 426 1.925222e-04
## GO:0060562 epithelial tube morphogenesis 0.3565320 257 5.932837e-04
## GO:0035295 tube development 0.3565320 391 5.953254e-04
##
## $less
## p.geomean stat.mean p.val
## GO:0048285 organelle fission 1.536227e-15 -8.063910 1.536227e-15
## GO:0000280 nuclear division 4.286961e-15 -7.939217 4.286961e-15
## GO:0007067 mitosis 4.286961e-15 -7.939217 4.286961e-15
12
## GO:0000087 M phase of mitotic cell cycle 1.169934e-14 -7.797496 1.169934e-14
## GO:0007059 chromosome segregation 2.028624e-11 -6.878340 2.028624e-11
## GO:0000236 mitotic prometaphase 1.729553e-10 -6.695966 1.729553e-10
## q.val set.size exp1
## GO:0048285 organelle fission 5.841698e-12 376 1.536227e-15
## GO:0000280 nuclear division 5.841698e-12 352 4.286961e-15
## GO:0007067 mitosis 5.841698e-12 352 4.286961e-15
## GO:0000087 M phase of mitotic cell cycle 1.195672e-11 362 1.169934e-14
## GO:0007059 chromosome segregation 1.658603e-08 142 2.028624e-11
## GO:0000236 mitotic prometaphase 1.178402e-07 84 1.729553e-10
##
## $stats
## stat.mean exp1
## GO:0007156 homophilic cell adhesion 3.824205 3.824205
## GO:0002009 morphogenesis of an epithelium 3.653886 3.653886
## GO:0048729 tissue morphogenesis 3.643242 3.643242
## GO:0007610 behavior 3.565432 3.565432
## GO:0060562 epithelial tube morphogenesis 3.261376 3.261376
## GO:0035295 tube development 3.253665 3.253665
Thegreaterresultsshowenrichedprocessesforup-regulatedgenes(e.g.epitheliumdevelopment, homophilic
cell adhesion), and thelessresults show enriched cell-cycle related processes (M phase, mitosis, etc.), which
matches the KEGG cell cycle signal.
Section 4 – Reactome Analysis
In this part I export a list of significant genes and upload it to theReactomewebsite to do over-
representation analysis.
# Genes significant at padj <= 0.05
sig_genes <- res[res$padj<=0.05& !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:",length(sig_genes)))
## [1] "Total number of significant genes: 8147"
write.table(sig_genes,
file = "significant_genes.txt",
row.names = FALSE,
col.names = FALSE,
quote = FALSE)
I then go tohttps://reactome.org/PathwayBrowser/#TOOL=AT, uploadsignificant_genes.txt, select
“Project to Humans”, and run the analysis.
Q8 – Reactome question:
“What pathway has the most significant ‘Entities p-value’? Do the most significant pathways
listed match your previous KEGG results? What factors could cause differences between the two
methods?”
13
My answer (Q8):
•In Reactome, the pathway with the most significantEntities p-valueis acell-cycle–related path-
way(for example, “Cell Cycle, Mitotic”), which lines up with what we saw from KEGG where the
Cell cycle (hsa04110)pathway was the top down-regulated pathway.
•So yes, the most significant Reactome pathways generallymatch the KEGG resultsin the sense
that both are strongly highlighting cell-cycle and DNA replication processes as being heavily perturbed
by HOXA1 knockdown.
•Differences between the methods can come from:
–Differentpathway databases and curation(KEGG vs Reactome have different gene sets and
boundaries for each pathway).
–Differentbackground gene universeor how each tool defines the set of “testable” genes.
–Slightly differentstatistical models, multiple-testing corrections, or enrichment metrics.
–Differences ingene ID mapping(Entrez vs Ensembl vs gene symbols) and which genes end up
being recognized by each tool.
Section 5 – GO Online Enrichment (Optional)
Here I take the same significant gene list and use the online GO enrichment portal to look specifically at
Biological Processterms.
Steps I follow:
1. Go to the GO enrichment website:
http://www.geneontology.org/page/go-enrichment-analysis
2. Paste thesignificant_genes.txtgene list.
3. Choose“Biological process”and“Homo sapiens”.
4. Submit and look at the top GO terms.
Q9 – GO online question:
“What pathway (GO term) has the most significant p-value? Do the most significant pathways
listed match your previous KEGG results? What factors could cause differences between the two
methods?”
My answer (Q9):
•ThemostsignificantGOBiologicalProcesstermsarealsocell-cycle–related, like“Mphase”, “mitotic
cell cycle”, or “nuclear division”. These GO terms match the same biological story as KEGG and
Reactome: HOXA1 knockdown strongly disrupts the cell cycle.
•So again, the top GO BP termsagree wellwith the KEGG and Reactome results: all three methods
say that genes involved in cell cycle and mitosis are heavily enriched among the differentially expressed
genes.
•Differences between GO and KEGG/Reactome can be due to:
–GO terms arehierarchical and more granular, so a single KEGG pathway may correspond
to many overlapping GO BP terms.
–Differentannotation coverage: some genes are well annotated in GO but not in KEGG (or vice
versa).
14
–Differentsize and structureof gene sets, which changes the enrichment statistics even with the
same gene list.
–SlightdifferencesinhowtheonlineGOtooldefinesthebackgroundgeneuniverseandperforms
multiple-testing correction.
Session Information (Reproducibility)
Finally, I include the session info so the analysis is reproducible.
sessionInfo()
## R version 4.5.1 (2025-06-13)
## Platform: aarch64-apple-darwin20
## Running under: macOS Ventura 13.3
##
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib; LAPACK version 3.12.1
##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## time zone: America/Los_Angeles
## tzcode source: internal
##
## attached base packages:
## [1] stats4 stats graphics grDevices utils datasets methods
## [8] base
##
## other attached packages:
## [1] gageData_2.48.0 gage_2.60.0
## [3] pathview_1.50.0 org.Hs.eg.db_3.22.0
## [5] AnnotationDbi_1.72.0 DESeq2_1.50.0
## [7] SummarizedExperiment_1.40.0 Biobase_2.70.0
## [9] MatrixGenerics_1.22.0 matrixStats_1.5.0
## [11] GenomicRanges_1.62.0 Seqinfo_1.0.0
## [13] IRanges_2.44.0 S4Vectors_0.48.0
## [15] BiocGenerics_0.56.0 generics_0.1.4
##
## loaded via a namespace (and not attached):
## [1] KEGGREST_1.50.0 gtable_0.3.6 xfun_0.52
## [4] ggplot2_4.0.0 lattice_0.22-7 bitops_1.0-9
## [7] vctrs_0.6.5 tools_4.5.1 parallel_4.5.1
## [10] tibble_3.3.0 RSQLite_2.4.3 blob_1.2.4
## [13] pkgconfig_2.0.3 Matrix_1.7-3 RColorBrewer_1.1-3
## [16] S7_0.2.0 graph_1.88.0 lifecycle_1.0.4
## [19] compiler_4.5.1 farver_2.1.2 Biostrings_2.78.0
## [22] codetools_0.2-20 htmltools_0.5.8.1 RCurl_1.98-1.17
## [25] yaml_2.3.10 GO.db_3.22.0 pillar_1.11.1
## [28] crayon_1.5.3 BiocParallel_1.44.0 DelayedArray_0.36.0
15
## [31] cachem_1.1.0 abind_1.4-8 tidyselect_1.2.1
## [34] locfit_1.5-9.12 digest_0.6.37 dplyr_1.1.4
## [37] fastmap_1.2.0 grid_4.5.1 cli_3.6.5
## [40] SparseArray_1.10.1 magrittr_2.0.3 S4Arrays_1.10.0
## [43] XML_3.99-0.19 scales_1.4.0 bit64_4.6.0-1
## [46] rmarkdown_2.29 XVector_0.50.0 httr_1.4.7
## [49] bit_4.6.0 png_0.1-8 memoise_2.0.1
## [52] evaluate_1.0.4 knitr_1.50 rlang_1.1.6
## [55] Rcpp_1.1.0 glue_1.8.0 DBI_1.2.3
## [58] Rgraphviz_2.54.0 KEGGgraph_1.70.0 rstudioapi_0.17.1
## [61] R6_2.6.1
16
