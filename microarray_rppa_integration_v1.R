# Microarray-RPPA Integration 
# 1/14/20

library(tidyverse)
library(data.table)
library(limma)

# microarray DE
microarray <- fread("~/microarray/MCF7_Y537S_expression_values.csv", header = T, stringsAsFactors = F)

## filter out columns don't want
### Only doing MCF7 parental (P), AR, IR and Y537s P, AR, IR
microarray <- microarray[,c(1,4,9:20)]

## Rename Rae columns with sample names
colnames(microarray)[3:14] <- c("P_Y537S_1", "P_Y537S_2","Y537S_AR_1", "Y537S_AR_2", "Y537S_IR_1","Y537S_IR_2",
                                "P_MCF7_1", "P_MCF7_2", "MCF7_AR_1","MCF7_AR_2", "MCF7_IR_1", "MCF7_IR_2") 

## Filter out probes that don't have an associated symbol
## Multiple probes can be for the same gene
microarray <- microarray %>% filter(Symbol != "")

## MCF7 expression data
mcf7 <- microarray %>% select(P_MCF7_1, P_MCF7_2, MCF7_AR_1, MCF7_AR_2, MCF7_IR_1, MCF7_IR_2)
rownames(mcf7) <- microarray$ProbeID

## calculate mean for each sample type for each gene
mcf7_mean <- apply(mcf7[,1:2], 1, mean)
mcf7_mean <- cbind(mcf7_mean, apply(mcf7[,3:4], 1, mean))
mcf7_mean <- as.data.frame(cbind(mcf7_mean, apply(mcf7[,5:6], 1, mean)))
colnames(mcf7_mean) <- c('control', 'AR', 'IR')
mcf7_mean$gene <- microarray$Symbol

## DE with limma
### P = 1; AR = 2; IR = 3
ContrastModel <- model.matrix(~0+factor(c(1,1,2,2,3,3)))
colnames(ContrastModel) <- c("P_MCF7", "AR_MCF7", "IR_MCF7")
fit_1 <- lmFit(mcf7, ContrastModel)
mcf7_contrast <- makeContrasts(AR_MCF7-P_MCF7,IR_MCF7-P_MCF7, AR_MCF7-IR_MCF7, levels = ContrastModel)
fit_2 <- contrasts.fit(fit_1, mcf7_contrast)
fit_2 <- eBayes(fit_2)

## AR vs. parental
AR_P_DE <- topTable(fit_2, coef=1, adjust="BH", nrow(mcf7))
AR_P_DE <- AR_P_DE[order(rownames(AR_P_DE)),]
AR_P_DE$gene <- microarray$Symbol

## IR vs. parental
IR_P_DE <- topTable(fit_2, coef=2, adjust="BH", nrow(mcf7))
IR_P_DE <- IR_P_DE[order(rownames(IR_P_DE)),]
IR_P_DE$gene <- microarray$Symbol

## AR vs. IR
AR_IR_DE <- topTable(fit_2, coef=3, adjust="BH", nrow(mcf7))
AR_IR_DE <- AR_IR_DE[order(rownames(AR_IR_DE)),]
AR_IR_DE$gene <- microarray$Symbol

'------------------------------------------------------------------------------------------------------------'
# RPPA DE
rppa <- fread("~/rppa/rppa/RPPA_data_analysis_v1_011620.csv", header = T, stringsAsFactors = F)
rppa <- as.data.frame(rppa[,1:20])

rppa_expression <- rppa[,c(5:7,9:11,13:15)]
rownames(rppa_expression) <- paste('ab', 1:nrow(rppa), sep = '_')
rppa_expression <- t(rppa_expression)

rppa_test <- rppa[,c(1:4,8,12,16)]
rppa_test$row_name <- paste('ab', 1:nrow(rppa), sep = '_')

### AR vs. P
rppa_test$AR_P_log2FC <- log2(rppa$mean_AR) - log2(rppa$mean_control)
rppa_test$AR_P_p_val <- sapply(1:ncol(rppa_expression), function(x) wilcox.test(rppa_expression[4:6,x], rppa_expression[1:3,x])$p.value)
rppa_test$AR_P_adj_p_val <- p.adjust(rppa_test$AR_P_p_val, method = 'BH')

### IR vs. P
rppa_test$IR_P_log2FC <- log2(rppa$mean_IR) - log2(rppa$mean_control)
rppa_test$IR_P_p_val <- sapply(1:ncol(rppa_expression), function(x) wilcox.test(rppa_expression[7:9,x], rppa_expression[1:3,x])$p.value)
rppa_test$IR_P_adj_p_val <- p.adjust(rppa_test$IR_P_p_val, method = 'BH')

'------------------------------------------------------------------------------------------------------------'
# Transformation of RPPA/microarray mean

## check for multiple expression values for one gene
sum(duplicated(rppa_test$gene))
sum(duplicated(microarray$Symbol))

## find common genes between rppa and microarray datasets
rppa_genes <- unique(rppa_test$gene)
microarray_genes <- unique(microarray$Symbol)
common_genes <- rppa_genes[rppa_genes %in% microarray_genes]

rppa_phosph_genes <- unique(rppa_test$gene[rppa_test$phosphorylated])
common_phosph_genes <- common_genes[common_genes %in% rppa_phosph_genes]

microarray_duplicate_genes <- unique(microarray$Symbol[duplicated(microarray$Symbol)])

## create rppa and microarray mean data.frame with common genes
rppa_micro <- data.frame()

for (current_gene in common_genes) {
  temp_micro <- mcf7_mean %>% filter(gene == current_gene) %>% select(c(control, AR, IR)) %>%
    rename(micro_control = control, micro_AR = AR, micro_IR = IR)
  temp_rppa <- rppa_test %>% filter(gene == current_gene) %>% select(c(mean_control, mean_AR, mean_IR, phosphorylated)) %>%
    rename(rppa_control = mean_control, rppa_AR = mean_AR, rppa_IR = mean_IR)
  
  temp_merge <- merge(temp_micro, temp_rppa)
  temp_merge$gene <- rep(current_gene, times = nrow(temp_merge))
  
  rppa_micro <- rbind(rppa_micro, temp_merge)
}

## transform (rppa to microarray)
### control
control_lm <- lm(data = rppa_micro, micro_control ~ rppa_control)
control_b <- summary(control_lm)$coefficients[1,1]
control_m <- summary(control_lm)$coefficients[2,1]
rppa_micro$rppa_trans_control <- rppa_micro$rppa_control*control_m + control_b
plot(data = rppa_micro, rppa_trans_control ~ rppa_control)
abline(b = control_m, a = control_b)

### AR
AR_lm <- lm(data = rppa_micro, micro_AR ~ rppa_AR)
AR_b <- summary(AR_lm)$coefficients[1,1]
AR_m <- summary(AR_lm)$coefficients[2,1]
rppa_micro$rppa_trans_AR <- rppa_micro$rppa_AR*AR_m + AR_b
plot(data = rppa_micro, rppa_trans_AR ~ rppa_AR)
abline(b = AR_m, a = AR_b)

### IR
IR_lm <- lm(data = rppa_micro, micro_IR ~ rppa_IR)
IR_b <- summary(IR_lm)$coefficients[1,1]
IR_m <- summary(IR_lm)$coefficients[2,1]
rppa_micro$rppa_trans_IR <- rppa_micro$rppa_IR*IR_m + IR_b
plot(data = rppa_micro, rppa_trans_IR ~ rppa_IR)
abline(b = IR_m, a = IR_b)

### nicer plot
ggplot(data = rppa_micro)+
  geom_point(aes(x = rppa_control, y = micro_control)) +
  geom_point(aes(x = rppa_control, y = rppa_trans_control), color = "blue") +
  geom_abline(slope = control_m, intercept = control_b, color = "red") +
  xlab("Control RPPA Mean") +
  ylab("") 

'----------------------------------------------------------------------------'
# Calculate log2FC for transformed RPPA