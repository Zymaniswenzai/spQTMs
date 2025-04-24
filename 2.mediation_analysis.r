library(mediation)

#Read gene list related to AS regulators in the GeneCards database
gene=read.csv('AS_regulator_list.csv')

# Define cancer types and working directories
cancer <- c('BLCA','BRCA','COAD','ESCA','HNSC', 
            'KIRC','KIRP','LIHC','LUAD','LUSC',
            'PAAD','SKCM','STAD','THCA','CESC',
            'LAML','LGG','PCPG','PRAD','SARC','UCEC')

wd <- c('~/eQTM/',
        '~/spQTM/')

# Prepare *_list.txt files
rec <- data.frame(cancer=NA, QTM=NA, gene=NA, splicing=NA, time=NA)

for(i in seq_along(cancer)){
  eqtm <- read.table(paste0(wd[1], cancer[i], '_CisSig.txt'), sep='\t', header=TRUE)
  sqtm <- read.table(paste0(wd[2], cancer[i], '_TransSig.txt'), sep='\t', header=TRUE)
  eqtm <- merge(eqtm, gene, by='gene')  # <— MODIFY THIS LINE based on your context

  eqtm <- eqtm[,1:2]
  eqtm1 <- unique(data.frame(snps=eqtm$snps))
  sqtm <- merge(eqtm1, sqtm, by='snps')
  tmp <- merge(eqtm, sqtm, by='snps')
  
  rec[i,] <- c(cancer[i],
               length(unique(tmp$snps)),
               length(unique(tmp$gene.x)),
               length(unique(tmp$gene.y)),
               nrow(tmp))
  
  tmp <- tmp[,1:3]
  colnames(tmp) <- c('cpg','gene','splicing')
  write.table(tmp, paste0(cancer[i], '_list.txt'), sep='\t', quote=FALSE, row.names=FALSE)
}

# Mediation analysis
setwd('~/05.mediation/')

col_names <- c(
  "cpg", "splicing", "gene", 
  "d_avg", "d_avg_p", "d_avg_ci_l", "d_avg_ci_u",
  "z_avg", "z_avg_p", "z_avg_ci_l", "z_avg_ci_u",
  "tau_coef", "tau_p", "tau_ci_l", "tau_ci_u",
  "n_avg", "n_avg_p", "n_avg_ci_l", "n_avg_ci_u"
)

for(i in seq_along(cancer)){
  list_data <- read.table(paste0(cancer[i], '_list.txt'), sep='\t', header=TRUE)
  data <- read.table(paste0(cancer[i], '_mediation.txt'), sep='\t', header=TRUE, check.names=FALSE)
  if ("Row.names" %in% colnames(data)) data$Row.names <- NULL

  res <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(res) <- col_names

  for(j in 1:nrow(list_data)){
    cpg <- list_data[j, "cpg"]
    gene <- list_data[j, "gene"]
    splicing <- list_data[j, "splicing"]

    selected_cols <- cpg %in% colnames(data) & gene %in% colnames(data) & splicing %in% colnames(data)
    if (!selected_cols) next

    tmp <- data[, c(cpg, gene, splicing)]
    colnames(tmp) <- c("Predictor", "Mediator", "Outcome")

    model_mediator <- lm(Mediator ~ Predictor, data = tmp)
    model_outcome  <- lm(Outcome ~ Predictor + Mediator, data = tmp)
    set.seed(123)
    result <- mediate(model_mediator, model_outcome, treat='Predictor', mediator='Mediator', bootstrap=1000)
    summary_result <- summary(result)

    row_result <- data.frame(
      cpg = cpg,
      splicing = splicing,
      gene = gene,
      d_avg = summary_result$d.avg,
      d_avg_p = summary_result$d.avg.p,
      d_avg_ci_l = summary_result$d.avg.ci[1],
      d_avg_ci_u = summary_result$d.avg.ci[2],
      z_avg = summary_result$z.avg,
      z_avg_p = summary_result$z.avg.p,
      z_avg_ci_l = summary_result$z.avg.ci[1],
      z_avg_ci_u = summary_result$z.avg.ci[2],
      tau_coef = summary_result$tau.coef,
      tau_p = summary_result$tau.p,
      tau_ci_l = summary_result$tau.ci[1],
      tau_ci_u = summary_result$tau.ci[2],
      n_avg = summary_result$n.avg,
      n_avg_p = summary_result$n.avg.p,
      n_avg_ci_l = summary_result$n.avg.ci[1],
      n_avg_ci_u = summary_result$n.avg.ci[2],
      stringsAsFactors = FALSE
    )

    # Filter: p < 0.05
    if (row_result$d_avg_p < 0.05) {
      res <- rbind(res, row_result)
    }
  }

  write.table(res, paste0(cancer[i], '_mediation_result.txt'), sep='\t', quote=FALSE, row.names=FALSE)
}
