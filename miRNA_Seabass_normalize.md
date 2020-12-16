# miRNA Sea Bass normalization

###We want our final data in Excel
library(readxl)
library(writexl)

###For analysis, we use the package DESeq2
library("DESeq2")

###I have received an excel file with raw reads mapped to miRNAs called "Dla_annotation_count"
original_data <- read_excel("C:\\r_data\\Dla_annotation_count.xlsx",range = cell_cols("A:X"))


row_name <- original_data[["Anno_idx"]]
data <- original_data[ -c(1:8) ]

cts <- as.matrix(data)

rownames(cts) <- row_name

coldata  <- matrix("treated", length(colnames(cts)), 1)
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("condition")

###We have to define which samples are treated and which are untreated
coldata["MLT4","condition"] <- "untreated"
coldata["MLT3","condition"] <- "untreated"
coldata["MLT2","condition"] <- "untreated"
coldata["MLT1","condition"] <- "untreated"
coldata["MHT12","condition"] <- "treated"
coldata["MHT11","condition"] <- "treated"
coldata["MHT10","condition"] <- "treated"
coldata["MHT9","condition"] <- "treated"
coldata["FLT8","condition"] <- "untreated"
coldata["FLT7","condition"] <- "untreated"
coldata["FLT6","condition"] <- "untreated"
coldata["FLT5","condition"] <- "untreated"
coldata["FHT16","condition"] <- "treated"
coldata["FHT15","condition"] <- "treated"
coldata["FHT14","condition"] <- "treated"
coldata["FHT13","condition"] <- "treated"

coldata <- as.data.frame(coldata)

coldata$condition <- factor(coldata$condition)

###We have to normalize the reads
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- estimateSizeFactors(dds); 
normalData <- counts(dds, normalized=TRUE)



normalData <- as.data.frame(normalData)
class(normalData)

colnames(normalData) <- paste(colnames(normalData), "_normal", sep="")

normalData["Anno_idx"] <-rownames(normalData)

total <-  merge(original_data, normalData,by="Anno_idx", sort= FALSE)

###Our output, the excel file with normalized reads
write_xlsx(total, "C:\\r_data\\Dla_annotation_count_res.xlsx")




