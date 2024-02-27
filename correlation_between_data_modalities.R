library(readxl)
library(tidyverse)
library(glue)
library(ggplot2)
library(broom)
library(writexl)
library(corrplot)
library(ComplexHeatmap)

correlation_analysis <- function(data_mat, outname, cluster_k = 20){
  
  cytokine_cor <- cor(data_mat, method="spearman", use = "pairwise.complete.obs")
  p.mat <- cor.mtest(data_mat, method="spearman", use = "pairwise.complete.obs")
  
  pval_out <- p.mat$p
  dimnames(pval_out) <- dimnames(cytokine_cor)
  
  FDR_out <- sapply(c(1:ncol(pval_out)), function(x) p.adjust(pval_out[,x], method="BH"))
  colnames(FDR_out) <- colnames(pval_out)
  
  write_xlsx(data.frame(feature = row.names(cytokine_cor), cytokine_cor), path=glue("correlation_spearman_{outname}.xlsx"))
  write_xlsx(data.frame(feature = row.names(pval_out), pval_out), path=glue("correlation_spearman_pvalue_{outname}.xlsx"))
  write_xlsx(data.frame(feature = row.names(FDR_out), FDR_out), path=glue("correlation_spearman_FDR_{outname}.xlsx"))
  
  FDR_out_plot <- FDR_out
  FDR_out_plot[pval_out > 0.05] <- 1
  
  grDevices::cairo_pdf(glue("correlation_plot_{outname}.pdf"), width=12, height=12)
  corrplot(cytokine_cor, type="full", order="hclust", tl.col = "black", tl.srt = 90,
           p.mat = FDR_out_plot, sig.level = 0.25, insig = "blank", diag=TRUE,
           title = outname,
           mar=c(0,0,4,0))
  corrRect.hclust(cytokine_cor, k = cluster_k, lwd=3)
  dev.off()
  
  grDevices::cairo_pdf(glue("correlation_plot_{outname}_all_value.pdf"), width=12, height=12)
  corrplot(cytokine_cor, type="full", order="hclust", tl.col = "black", tl.srt = 90,
           p.mat = FDR_out_plot, sig.level = 1, insig = "blank", diag=TRUE,
           title = outname,
           mar=c(0,0,4,0))
  dev.off()
  
  grDevices::cairo_pdf(glue("correlation_plot_{outname}_hclust.pdf"), width=12, height=6)
  plot(hclust(as.dist(1-cytokine_cor)))
  dev.off()
  
  return(list(cytokine_cor, pval_out, FDR_out))
}

regression_analysis <- function(data_x, data_y, outname, FDR_cutoff = 0.25, pval_cutoff = 0.05){
  
  ##get t statistics
  stat_list <- lapply(colnames(data_y), function(d){
    data_y_select <- data_y[,d]
    
    stat_lm <- lapply(colnames(data_x), function(x){
      data_x_select <- data_x[,x]
      data_fit <- data.frame(y=data_y_select, x=data_x_select)
      fit_lm <- lm(y ~ x, data_fit)
      stat_lm <- tidy(fit_lm)
      return(data.frame(t_lm = stat_lm$statistic[2], p_lm = stat_lm$p.value[2]))
    })
    
    stat_lm <- Reduce(rbind, stat_lm)
    stat_lm$p_adj_lm <- p.adjust(stat_lm$p_lm, method="BH")
    row.names(stat_lm) <- colnames(data_x)
    return(stat_lm)
  })
  
  t_mat <- sapply(stat_list, function(x) x$t_lm)
  row.names(t_mat) <- colnames(data_x)
  colnames(t_mat) <- colnames(data_y)
  
  write_xlsx(data.frame(cytokine=row.names(t_mat), as.data.frame(t_mat)), glue("t_lm_{outname}.xlsx"))
  
  p_mat <- sapply(stat_list, function(x) x$p_lm)
  row.names(p_mat) <- colnames(data_x)
  colnames(p_mat) <- colnames(data_y)
  
  write_xlsx(data.frame(cytokine=row.names(p_mat), as.data.frame(p_mat)), glue("pvalue_lm_{outname}.xlsx"))
  
  p_adj_mat <- sapply(stat_list, function(x) x$p_adj_lm)
  row.names(p_adj_mat) <- colnames(data_x)
  colnames(p_adj_mat) <- colnames(data_y)
  
  write_xlsx(data.frame(cytokine=row.names(p_adj_mat), as.data.frame(p_adj_mat)), glue("FDR_lm_{outname}.xlsx"))
  
  ##set cutoff for significant marking, require p-value <= pval_cutoff and FDR <= FDR_cutoff 
  p_adj_mat_mark <- p_adj_mat
  p_adj_mat_mark[p_adj_mat <= FDR_cutoff] <- "*"
  p_adj_mat_mark[p_adj_mat > FDR_cutoff] <- ""
  p_adj_mat_mark[p_mat > pval_cutoff] <- ""
  
  ht3 <- Heatmap(t_mat,
                 name = "t-statistic",
                 column_title = outname,
                 cell_fun = function(j, i, x, y, w, h, col) {
                   grid.text(p_adj_mat_mark[i,j], x, y)
                 })
  
  grDevices::cairo_pdf(glue("lm_t_heatmap_{outname}.pdf"), width=10, height=10)
  draw(ht3)
  dev.off()
  
  return(list(t_mat = t_mat, p_mat = p_mat, p_adj_mat = p_adj_mat))
}

##main
data_in <- read_xlsx("data_file.xlsx")
                      
##correlation plot SOT
data_SOT <- data_in %>%
  filter(`Disease/group type` == "Transplant COVID-19") %>%
  rename(cfmtDNA = `cfmtDNA (Copies/mL)`, 
         ncfDNA = `ncfDNA (Copies/mL)`
         )

data_mat_SOT <- log2(data_SOT[,-c(1:8)] + 1)

correlation_analysis_SOT <- correlation_analysis(data_mat_SOT, "SOT", 20)

##correlation plot Non-SOT
data_Non_SOT <- data_in %>%
  filter(`Disease/group type` == "COVID-19") %>%
  rename(cfmtDNA = `cfmtDNA (Copies/mL)`, 
         ncfDNA = `ncfDNA (Copies/mL)`
  )

data_mat_Non_SOT <- log2(data_Non_SOT[,-c(1:8)] + 1)

correlation_analysis_Non_SOT <- correlation_analysis(data_mat_Non_SOT, "Non_SOT", 20)

##regression analysis between cytokine and cfDNA for SOT patients using cytokine as x and cfDNA as y
regression_analysis_SOT <- regression_analysis(data_x = data_mat_SOT[,-c(1:25)], data_y = data_mat_SOT[,c(1:25)], "SOT")

##regression analysis between cytokine and cfDNA for Non_SOT patients using cytokine as x and cfDNA as y
regression_analysis_Non_SOT <- regression_analysis(data_x = data_mat_Non_SOT[,-c(1:25)], data_y = data_mat_Non_SOT[,c(1:25)], "Non_SOT")
