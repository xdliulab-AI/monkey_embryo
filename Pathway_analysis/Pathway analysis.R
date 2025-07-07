#Pathway analysis
library(ggplot2)
library(Matrix)
library(Seurat)
library(data.table)

obj <- readRDS("Gene-corr-AP/D03259G613_bin50_Ddistance.rds")

activin <- list("activin" = unique(c('INHBA','ACVR1B','ACVR2A','INHBA','ACVR1B','ACVR2B','INHBB','ACVR1B','ACVR2A','INHBB','ACVR1B','ACVR2B','INHBB','ACVR1C','ACVR2A','INHBB','ACVR1C','ACVR2B','INHBC','ACVR1B','ACVR2A','INHBC','ACVR1B','ACVR2B','INHBC','ACVR1C','ACVR2A','INHBC','ACVR1C','ACVR2B','INHBE','ACVR1B','ACVR2A','INHBE','ACVR1B','ACVR2B','INHBE','ACVR1C','ACVR2A','INHBE','ACVR1C','ACVR2B','INHBABB','ACVR1B','ACVR2A','INHBABB','ACVR1B','ACVR2B','INHBABB','ACVR1C','ACVR2A','INHB','ACVR1C','ACVR2B','INHABA','ACVR2A','','INHABB','ACVR2B')))
ptn <- list("PTN" = unique(c('PTN','PTPRZ1','PTN','SDC1','PTN','SDC2','PTN','SDC3','PTN','SDC4','PTN','ITGAV','ITGB3','PTN','NCL','PTN','ALK')))
spp1 <- list("SPP1" = unique(c('SPP1','CD44','SPP1','ITGAV','ITGB3','SPP1','ITGAV','ITGB1','SPP1','ITGAV','ITGB5','SPP1','ITGAV','ITGB6','SPP1','ITGA4','ITGB1','SPP1','ITGA9','ITGB1','SPP1','ITGA8','ITGB1','SPP1','ITGA5','ITGB1')))
ncwnt <- list("ncWNT" = unique(c('WNT5A','FZD1','WNT5A','FZD10','WNT5A','FZD2','WNT5A','FZD3','WNT5A','FZD4','WNT5A','FZD5','WNT5A','FZD6','WNT5A','FZD7','WNT5A','FZD8','WNT5A','FZD9','WNT5A','MCAM','WNT5B','FZD1','WNT5B','FZD10','WNT5B','FZD2','WNT5B','FZD3','WNT5B','FZD4','WNT5B','FZD5','WNT5B','FZD6','WNT5B','FZD7','WNT5B','FZD8','WNT5B','FZD9','WNT11','FZD1','WNT11','FZD10','WNT11','FZD2','WNT11','FZD3','WNT11','FZD4','WNT11','FZD5','WNT11','FZD6','WNT11','FZD7','WNT11','FZD8','WNT11','FZD9')))
sema3 <- list("SEMA3" = unique(c('SEMA3A','NRP1','PLXNA1','SEMA3A','NRP1','PLXNA2','SEMA3A','NRP1','PLXNA3','SEMA3A','NRP1','PLXNA4','SEMA3B','NRP1','PLXNA1','SEMA3B','NRP1','PLXNA2','SEMA3B','NRP1','PLXNA3','SEMA3B','NRP1','PLXNA4','SEMA3C','NRP1','PLXNA1','SEMA3C','NRP1','PLXNA2','SEMA3C','NRP1','PLXNA3','SEMA3C','NRP1','PLXNA4','SEMA3D','NRP1','PLXNA1','SEMA3D','NRP1','PLXNA2','SEMA3D','NRP1','PLXNA3','SEMA3D','NRP1','PLXNA4','SEMA3B','NRP2','PLXNA1','SEMA3B','NRP2','PLXNA2','SEMA3B','NRP2','PLXNA3','SEMA3B','NRP2','PLXNA4','SEMA3C','NRP2','PLXNA1','SEMA3C','NRP2','PLXNA2','SEMA3C','NRP2','PLXNA3','SEMA3C','NRP2','PLXNA4','SEMA3D','NRP2','PLXNA1','SEMA3D','NRP2','PLXNA2','SEMA3D','NRP2','PLXNA3','SEMA3D','NRP2','PLXNA4','SEMA3F','NRP2','PLXNA1','SEMA3F','NRP2','PLXNA2','SEMA3F','NRP2','PLXNA3','SEMA3F','NRP2','PLXNA4','SEMA3G','NRP2','PLXNA1','SEMA3G','NRP2','PLXNA2','SEMA3G','NRP2','PLXNA3','SEMA3G','NRP2','PLXNA4','SEMA3C','NRP1','NRP2','SEMA3C','PLXND1','SEMA3E','PLXND1')))
bmp <- list("BMP" = unique(c('BMP2','BMPR1A','ACVR2A','BMP2','BMPR1A','ACVR2B','BMP2','BMPR1A','BMPR2','BMP2','BMPR1B','ACVR2A','BMP2','BMPR1B','ACVR2B','BMP2','BMPR1B','BMPR2','BMP4','BMPR1A','ACVR2A','BMP4','BMPR1A','ACVR2B','BMP4','BMPR1A','BMPR2','BMP4','BMPR1B','ACVR2A','BMP4','BMPR1B','ACVR2B','BMP4','BMPR1B','BMPR2','GDF5','BMPR1A','ACVR2A','GDF5','BMPR1A','ACVR2B','GDF5','BMPR1A','BMPR2','GDF5','BMPR1B','ACVR2A','GDF5','BMPR1B','ACVR2B','GDF5','BMPR1B','BMPR2','GDF6','BMPR1A','ACVR2A','GDF6','BMPR1A','ACVR2B','GDF6','BMPR1A','BMPR2','GDF6','BMPR1B','ACVR2A','GDF6','BMPR1B','ACVR2B','GDF6','BMPR1B','BMPR2','GDF7','BMPR1A','ACVR2A','GDF7','BMPR1A','ACVR2B','GDF7','BMPR1A','BMPR2','GDF7','BMPR1B','ACVR2A','GDF7','BMPR1B','ACVR2B','GDF7','BMPR1B','BMPR2','BMP15','BMPR1B','BMPR2','BMP5','ACVR1','ACVR2A','BMP5','ACVR1','ACVR2B','BMP5','ACVR1','BMPR2','BMP5','BMPR1A','ACVR2A','BMP5','BMPR1A','ACVR2B','BMP5','BMPR1A','BMPR2','BMP5','BMPR1B','ACVR2A','BMP5','BMPR1B','ACVR2B','BMP5','BMPR1B','BMPR2','BMP6','ACVR1','ACVR2A','BMP6','ACVR1','ACVR2B','BMP6','ACVR1','BMPR2','BMP6','BMPR1A','ACVR2A','BMP6','BMPR1A','ACVR2B','BMP6','BMPR1A','BMPR2','BMP6','BMPR1B','ACVR2A','BMP6','BMPR1B','ACVR2B','BMP6','BMPR1B','BMPR2','BMP7','ACVR1','ACVR2A','BMP7','ACVR1','ACVR2B','BMP7','ACVR1','BMPR2','BMP7','BMPR1A','ACVR2A','BMP7','BMPR1A','ACVR2B','BMP7','BMPR1A','BMPR2','BMP7','BMPR1B','ACVR2A','BMP7','BMPR1B','ACVR2B','BMP7','BMPR1B','BMPR2','BMP8A','ACVR1','ACVR2A','BMP8A','ACVR1','ACVR2B','BMP8A','ACVR1','BMPR2','BMP8A','BMPR1A','ACVR2A','BMP8A','BMPR1A','ACVR2B','BMP8A','BMPR1A','BMPR2','BMP8A','BMPR1B','ACVR2A','BMP8A','BMPR1B','ACVR2B','BMP8A','BMPR1B','BMPR2','BMP8B','ACVR1','ACVR2A','BMP8B','ACVR1','ACVR2B','BMP8B','ACVR1','BMPR2','BMP8B','BMPR1A','ACVR2A','BMP8B','BMPR1A','ACVR2B','BMP8B','BMPR1A','BMPR2','BMP8B','BMPR1B','ACVR2A','BMP8B','BMPR1B','ACVR2B','BMP8B','BMPR1B','BMPR2')))
fgf <- list("FGF" = unique(c('FGF1','FGFR1','FGF1','FGFR2','FGF1','FGFR3','FGF1','FGFR4','FGF2','FGFR1','FGF2','FGFR2','FGF2','FGFR3','FGF2','FGFR4','FGF4','FGFR1','FGF4','FGFR2','FGF4','FGFR3','FGF4','FGFR4','FGF5','FGFR1','FGF5','FGFR2','FGF5','FGFR3','FGF5','FGFR4','FGF6','FGFR1','FGF6','FGFR2','FGF6','FGFR3','FGF6','FGFR4','FGF7','FGFR1','FGF7','FGFR2','FGF3','FGFR1','FGF3','FGFR2','FGF10','FGFR1','FGF10','FGFR2','FGF22','FGFR2','FGF8','FGFR1','FGF8','FGFR2','FGF8','FGFR3','FGF8','FGFR4','FGF17','FGFR1','FGF17','FGFR2','FGF17','FGFR3','FGF17','FGFR4','FGF18','FGFR1','FGF18','FGFR2','FGF18','FGFR3','FGF18','FGFR4','FGF9','FGFR1','FGF9','FGFR2','FGF9','FGFR3','FGF9','FGFR4','FGF20','FGFR1','FGF20','FGFR2','FGF20','FGFR3','FGF20','FGFR4','FGF16','FGFR1','FGF16','FGFR2','FGF16','FGFR3','FGF16','FGFR4','FGF15','FGFR1','FGF15','FGFR2','FGF15','FGFR3','FGF15','FGFR4','FGF21','FGFR1','FGF21','FGFR3','FGF23','FGFR1','FGF23','FGFR3','FGF23','FGFR4')))
IGF <- list("IGF" = unique(c('IGF1','IGF1R','IGF2','IGF1R','IGF2','IGF2R','IGF1','ITGAV','ITGB3','IGF1','ITGA6','ITGB4','IGF2','ITGAV','ITGB3','IGF2','ITGA6','ITGB4','IGFL3','IGFLR1')))
HH <- list("HH" = unique(c('DHH','PTCH1','IHH','PTCH1','SHH','PTCH1','DHH','PTCH2','IHH','PTCH2','SHH','PTCH2')))
cxcl <- list("CXCL" = unique(c('CXCL1','ACKR1','CXCL2','ACKR1','CXCL3','ACKR1','CXCL5','ACKR1','CXCL6','ACKR1','CXCL8','ACKR1','CXCL9','ACKR1','CXCL10','ACKR1','CXCL11','ACKR1','CXCL13','ACKR1','CXCL1','CXCR1','CXCL2','CXCR1','CXCL3','CXCR1','CXCL5','CXCR1','CXCL6','CXCR1','CXCL8','CXCR1','CXCL1','CXCR2','CXCL2','CXCR2','CXCL3','CXCR2','CXCL5','CXCR2','CXCL6','CXCR2','CXCL8','CXCR2','PPBP','CXCR2','PF4V1','CXCR3','CXCL9','CXCR3','CXCL10','CXCR3','CXCL11','CXCR3','CXCL13','CXCR3','PF4','CXCR3','CXCL12','CXCR4','CXCL11','ACKR3','CXCL12','ACKR3','CXCL13','CXCR5','CXCL16','CXCR6')))
grn <- list("GRN" = unique(c('GRN','SORT1')))
GDF <- list("GDF" = unique(c('GDF9','ACVR1B','BMPR2','GDF10','ACVR1B','ACVR2A','GDF11','TGFBR1','ACVR2A','GDF11','TGFBR1','ACVR2B','GDF11','ACVR1B','ACVR2B','GDF15','TGFBR2')))
NODAL <- list("NODAL" = unique(c('NODAL','ACVR1B','ACVR2A','NODAL','ACVR1C','ACVR2A','NODAL','ACVR1B','ACVR2B','NODAL','ACVR1C','ACVR2B','GDF1','ACVR1B','ACVR2A','GDF1','ACVR1B','ACVR2B','GDF1','ACVR1C','ACVR2A','GDF1','ACVR1C','ACVR2B','GDF3','ACVR1B','ACVR2A','GDF3','ACVR1B','ACVR2B','GDF3','ACVR1C','ACVR2A','GDF3','ACVR1C','ACVR2B')))
WNT <- list("WNT" = unique(c('WNT1','FZD1','LRP5','WNT1','FZD10','LRP5','WNT1','FZD2','LRP5','WNT1','FZD3','LRP5','WNT1','FZD4','LRP5','WNT1','FZD5','LRP5','WNT1','FZD6','LRP5','WNT1','FZD7','LRP5','WNT1','FZD8','LRP5','WNT1','FZD9','LRP5','WNT10A','FZD1','LRP5','WNT10A','FZD10','LRP5','WNT10A','FZD2','LRP5','WNT10A','FZD3','LRP5','WNT10A','FZD4','LRP5','WNT10A','FZD5','LRP5','WNT10A','FZD6','LRP5','WNT10A','FZD7','LRP5','WNT10A','FZD8','LRP5','WNT10A','FZD9','LRP5','WNT10B','FZD1','LRP5','WNT10B','FZD10','LRP5','WNT10B','FZD2','LRP5','WNT10B','FZD3','LRP5','WNT10B','FZD4','LRP5','WNT10B','FZD5','LRP5','WNT10B','FZD6','LRP5','WNT10B','FZD7','LRP5','WNT10B','FZD8','LRP5','WNT10B','FZD9','LRP5','WNT16','FZD1','LRP5','WNT16','FZD10','LRP5','WNT16','FZD2','LRP5','WNT16','FZD3','LRP5','WNT16','FZD4','LRP5','WNT16','FZD5','LRP5','WNT16','FZD6','LRP5','WNT16','FZD7','LRP5','WNT16','FZD8','LRP5','WNT16','FZD9','LRP5','WNT2','FZD1','LRP5','WNT2','FZD10','LRP5','WNT2','FZD2','LRP5','WNT2','FZD3','LRP5','WNT2','FZD4','LRP5','WNT2','FZD5','LRP5','WNT2','FZD6','LRP5','WNT2','FZD7','LRP5','WNT2','FZD8','LRP5','WNT2','FZD9','LRP5','WNT2B','FZD1','LRP5','WNT2B','FZD10','LRP5','WNT2B','FZD2','LRP5','WNT2B','FZD3','LRP5','WNT2B','FZD4','LRP5','WNT2B','FZD5','LRP5','WNT2B','FZD6','LRP5','WNT2B','FZD7','LRP5','WNT2B','FZD8','LRP5','WNT2B','FZD9','LRP5','WNT3','FZD1','LRP5','WNT3','FZD10','LRP5','WNT3','FZD2','LRP5','WNT3','FZD3','LRP5','WNT3','FZD4','LRP5','WNT3','FZD5','LRP5','WNT3','FZD6','LRP5','WNT3','FZD7','LRP5','WNT3','FZD8','LRP5','WNT3','FZD9','LRP5','WNT3A','FZD1','LRP5','WNT3A','FZD10','LRP5','WNT3A','FZD2','LRP5','WNT3A','FZD3','LRP5','WNT3A','FZD4','LRP5','WNT3A','FZD5','LRP5','WNT3A','FZD6','LRP5','WNT3A','FZD7','LRP5','WNT3A','FZD8','LRP5','WNT3A','FZD9','LRP5','WNT4','FZD1','LRP5','WNT4','FZD10','LRP5','WNT4','FZD2','LRP5','WNT4','FZD3','LRP5','WNT4','FZD4','LRP5','WNT4','FZD5','LRP5','WNT4','FZD6','LRP5','WNT4','FZD7','LRP5','WNT4','FZD8','LRP5','WNT4','FZD9','LRP5','WNT6','FZD1','LRP5','WNT6','FZD10','LRP5','WNT6','FZD2','LRP5','WNT6','FZD3','LRP5','WNT6','FZD4','LRP5','WNT6','FZD5','LRP5','WNT6','FZD6','LRP5','WNT6','FZD7','LRP5','WNT6','FZD8','LRP5','WNT6','FZD9','LRP5','WNT7A','FZD1','LRP5','WNT7A','FZD10','LRP5','WNT7A','FZD2','LRP5','WNT7A','FZD3','LRP5','WNT7A','FZD4','LRP5','WNT7A','FZD5','LRP5','WNT7A','FZD6','LRP5','WNT7A','FZD7','LRP5','WNT7A','FZD8','LRP5','WNT7A','FZD9','LRP5','WNT7B','FZD1','LRP5','WNT7B','FZD10','LRP5','WNT7B','FZD2','LRP5','WNT7B','FZD3','LRP5','WNT7B','FZD4','LRP5','WNT7B','FZD5','LRP5','WNT7B','FZD6','LRP5','WNT7B','FZD7','LRP5','WNT7B','FZD8','LRP5','WNT7B','FZD9','LRP5','WNT8A','FZD1','LRP5','WNT8A','FZD10','LRP5','WNT8A','FZD2','LRP5','WNT8A','FZD3','LRP5','WNT8A','FZD4','LRP5','WNT8A','FZD5','LRP5','WNT8A','FZD6','LRP5','WNT8A','FZD7','LRP5','WNT8A','FZD8','LRP5','WNT8A','FZD9','LRP5','WNT8B','FZD1','LRP5','WNT8B','FZD10','LRP5','WNT8B','FZD2','LRP5','WNT8B','FZD3','LRP5','WNT8B','FZD4','LRP5','WNT8B','FZD5','LRP5','WNT8B','FZD6','LRP5','WNT8B','FZD7','LRP5','WNT8B','FZD8','LRP5','WNT8B','FZD9','LRP5','WNT9A','FZD1','LRP5','WNT9A','FZD10','LRP5','WNT9A','FZD2','LRP5','WNT9A','FZD3','LRP5','WNT9A','FZD4','LRP5','WNT9A','FZD5','LRP5','WNT9A','FZD6','LRP5','WNT9A','FZD7','LRP5','WNT9A','FZD8','LRP5','WNT9A','FZD9','LRP5','WNT9B','FZD1','LRP5','WNT9B','FZD10','LRP5','WNT9B','FZD2','LRP5','WNT9B','FZD3','LRP5','WNT9B','FZD4','LRP5','WNT9B','FZD5','LRP5','WNT9B','FZD6','LRP5','WNT9B','FZD7','LRP5','WNT9B','FZD8','LRP5','WNT9B','FZD9','LRP5','WNT1','FZD1','LRP6','WNT1','FZD10','LRP6','WNT1','FZD2','LRP6','WNT1','FZD3','LRP6','WNT1','FZD4','LRP6','WNT1','FZD5','LRP6','WNT1','FZD6','LRP6','WNT1','FZD7','LRP6','WNT1','FZD8','LRP6','WNT1','FZD9','LRP6','WNT10A','FZD1','LRP6','WNT10A','FZD10','LRP6','WNT10A','FZD2','LRP6','WNT10A','FZD3','LRP6','WNT10A','FZD4','LRP6','WNT10A','FZD5','LRP6','WNT10A','FZD6','LRP6','WNT10A','FZD7','LRP6','WNT10A','FZD8','LRP6','WNT10A','FZD9','LRP6','WNT10B','FZD1','LRP6','WNT10B','FZD10','LRP6','WNT10B','FZD2','LRP6','WNT10B','FZD3','LRP6','WNT10B','FZD4','LRP6','WNT10B','FZD5','LRP6','WNT10B','FZD6','LRP6','WNT10B','FZD7','LRP6','WNT10B','FZD8','LRP6','WNT10B','FZD9','LRP6','WNT16','FZD1','LRP6','WNT16','FZD10','LRP6','WNT16','FZD2','LRP6','WNT16','FZD3','LRP6','WNT16','FZD4','LRP6','WNT16','FZD5','LRP6','WNT16','FZD6','LRP6','WNT16','FZD7','LRP6','WNT16','FZD8','LRP6','WNT16','FZD9','LRP6','WNT2','FZD1','LRP6','WNT2','FZD10','LRP6','WNT2','FZD2','LRP6','WNT2','FZD3','LRP6','WNT2','FZD4','LRP6','WNT2','FZD5','LRP6','WNT2','FZD6','LRP6','WNT2','FZD7','LRP6','WNT2','FZD8','LRP6','WNT2','FZD9','LRP6','WNT2B','FZD1','LRP6','WNT2B','FZD10','LRP6','WNT2B','FZD2','LRP6','WNT2B','FZD3','LRP6','WNT2B','FZD4','LRP6','WNT2B','FZD5','LRP6','WNT2B','FZD6','LRP6','WNT2B','FZD7','LRP6','WNT2B','FZD8','LRP6','WNT2B','FZD9','LRP6','WNT3','FZD1','LRP6','WNT3','FZD10','LRP6','WNT3','FZD2','LRP6','WNT3','FZD3','LRP6','WNT3','FZD4','LRP6','WNT3','FZD5','LRP6','WNT3','FZD6','LRP6','WNT3','FZD7','LRP6','WNT3','FZD8','LRP6','WNT3','FZD9','LRP6','WNT3A','FZD1','LRP6','WNT3A','FZD10','LRP6','WNT3A','FZD2','LRP6','WNT3A','FZD3','LRP6','WNT3A','FZD4','LRP6','WNT3A','FZD5','LRP6','WNT3A','FZD6','LRP6','WNT3A','FZD7','LRP6','WNT3A','FZD8','LRP6','WNT3A','FZD9','LRP6','WNT4','FZD1','LRP6','WNT4','FZD10','LRP6','WNT4','FZD2','LRP6','WNT4','FZD3','LRP6','WNT4','FZD4','LRP6','WNT4','FZD5','LRP6','WNT4','FZD6','LRP6','WNT4','FZD7','LRP6','WNT4','FZD8','LRP6','WNT4','FZD9','LRP6','WNT6','FZD1','LRP6','WNT6','FZD10','LRP6','WNT6','FZD2','LRP6','WNT6','FZD3','LRP6','WNT6','FZD4','LRP6','WNT6','FZD5','LRP6','WNT6','FZD6','LRP6','WNT6','FZD7','LRP6','WNT6','FZD8','LRP6','WNT6','FZD9','LRP6','WNT7A','FZD1','LRP6','WNT7A','FZD10','LRP6','WNT7A','FZD2','LRP6','WNT7A','FZD3','LRP6','WNT7A','FZD4','LRP6','WNT7A','FZD5','LRP6','WNT7A','FZD6','LRP6','WNT7A','FZD7','LRP6','WNT7A','FZD8','LRP6','WNT7A','FZD9','LRP6','WNT7B','FZD1','LRP6','WNT7B','FZD10','LRP6','WNT7B','FZD2','LRP6','WNT7B','FZD3','LRP6','WNT7B','FZD4','LRP6','WNT7B','FZD5','LRP6','WNT7B','FZD6','LRP6','WNT7B','FZD7','LRP6','WNT7B','FZD8','LRP6','WNT7B','FZD9','LRP6','WNT8A','FZD1','LRP6','WNT8A','FZD10','LRP6','WNT8A','FZD2','LRP6','WNT8A','FZD3','LRP6','WNT8A','FZD4','LRP6','WNT8A','FZD5','LRP6','WNT8A','FZD6','LRP6','WNT8A','FZD7','LRP6','WNT8A','FZD8','LRP6','WNT8A','FZD9','LRP6','WNT8B','FZD1','LRP6','WNT8B','FZD10','LRP6','WNT8B','FZD2','LRP6','WNT8B','FZD3','LRP6','WNT8B','FZD4','LRP6','WNT8B','FZD5','LRP6','WNT8B','FZD6','LRP6','WNT8B','FZD7','LRP6','WNT8B','FZD8','LRP6','WNT8B','FZD9','LRP6','WNT9A','FZD1','LRP6','WNT9A','FZD10','LRP6','WNT9A','FZD2','LRP6','WNT9A','FZD3','LRP6','WNT9A','FZD4','LRP6','WNT9A','FZD5','LRP6','WNT9A','FZD6','LRP6','WNT9A','FZD7','LRP6','WNT9A','FZD8','LRP6','WNT9A','FZD9','LRP6','WNT9B','FZD1','LRP6','WNT9B','FZD10','LRP6','WNT9B','FZD2','LRP6','WNT9B','FZD3','LRP6','WNT9B','FZD4','LRP6','WNT9B','FZD5','LRP6','WNT9B','FZD6','LRP6','WNT9B','FZD7','LRP6','WNT9B','FZD8','LRP6','WNT9B','FZD9','LRP6')))

obj <- AddModuleScore(obj,
                      features = activin,
                      ctrl = 100,
                      name = "Activin")
obj <- AddModuleScore(obj,
                      features = ptn,
                      ctrl = 100,
                      name = "PTN")
obj <- AddModuleScore(obj,
                      features = spp1,
                      ctrl = 100,
                      name = "SPP1")
obj <- AddModuleScore(obj,
                      features = ncwnt,
                      ctrl = 100,
                      name = "ncWNT")
obj <- AddModuleScore(obj,
                      features = sema3,
                      ctrl = 100,
                      name = "SEMA3")
obj <- AddModuleScore(obj,
                      features = bmp,
                      ctrl = 100,
                      name = "BMP")
obj <- AddModuleScore(obj,
                      features = fgf,
                      ctrl = 100,
                      name = "FGF")
obj <- AddModuleScore(obj,
                      features = IGF,
                      ctrl = 100,
                      name = "IGF")
obj <- AddModuleScore(obj,
                      features = HH,
                      ctrl = 100,
                      name = "HH")
obj <- AddModuleScore(obj,
                      features = cxcl,
                      ctrl = 100,
                      name = "CXCL")
obj <- AddModuleScore(obj,
                      features = grn,
                      ctrl = 100,
                      name = "GRN")
obj <- AddModuleScore(obj,
                      features = GDF,
                      ctrl = 100,
                      name = "GDF")
obj <- AddModuleScore(obj,
                      features = NODAL,
                      ctrl = 100,
                      name = "NODAL")
obj <- AddModuleScore(obj,
                      features = WNT,
                      ctrl = 100,
                      name = "WNT")


# 安装并加载dplyr包（如果未安装）
library(dplyr)

# 使用dplyr包的group_by和summarize函数按group列对value1、value2、value3列求平均值
averages_dplyr <- modulescore %>%
  group_by(Aportion_id) %>%
  summarize(across(c(Activin1,PTN1,SPP11,ncWNT1,SEMA31,BMP1,FGF1,IGF1,HH1,CXCL1,GRN1,GDF1,NODAL1,WNT1), mean))

averages_dplyr <- as.data.frame(averages_dplyr)
rownames(averages_dplyr) <- averages_dplyr$Aportion_id
averages_dplyrsub <- averages_dplyr[,-1]
averages_dplyrsubt <- as.data.frame(t(averages_dplyrsub))

#'Activin1','PTN1','SPP11','ncWNT1','SEMA31','BMP1','FGF1','IGF1','HH1','CXCL1','GRN1','GDF1','NODAL1','WNT1'
ggplot(averages_dplyr,aes(Aportion_id,Activin1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                                  color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,PTN1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,SPP11))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                               color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,ncWNT1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                                color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,SEMA31))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                                color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,BMP1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,FGF1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,IGF1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))


ggplot(averages_dplyr,aes(Aportion_id,HH1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                             color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,CXCL1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                               color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,GRN1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,GDF1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,NODAL1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                                color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))

ggplot(averages_dplyr,aes(Aportion_id,WNT1))+geom_point(size=0.5)+geom_smooth(method = 'loess',se=TRUE,
                                                                              color="black", formula = y ~ x)+theme_bw()+scale_fill_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+
  scale_color_manual(values = c("#5c6672", "#35aae0","#eb5b26"))+theme(axis.text = element_text(color = "black"))