rm(list = ls())
library(maftools)
options(stringsAsFactors = F)
setwd("/home/eunji/proj/wgstest/")

#annovar
annovar.maf <- annovarToMaf(annovar = "./project/7.annotation/annovar/annovar_merge.vcf",
                        refBuild = 'hg38',
                        tsbCol = 'Tumor_Sample_Barcode',
                        table = 'refGene',
                        MAFobj = T)
table(annovar.maf@data[annovar.maf@data$Tumor_Sample_Barcode=='15190-MLS-D',]$Variant_Type)

maf=annovar.maf
unique(maf@data$Tumor_Sample_Barcode)
getSampleSummary(maf)
getGeneSummary(maf)
getFields(maf)

plotmafSummary(maf = maf,
               rmOutlier = T,
               showBarcodes = T,
               textSize = 0.4,
               addStat = 'median',
               dashboard = T,
               titvRaw = F)
