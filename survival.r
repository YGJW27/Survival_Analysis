# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#The_lung_dataset
library(survival)
library(survminer)

# load file
file_name <- "D:/luad_tcga_pan_can_atlas_2018_clinical_data.tsv"
rppa_name <- "D:/PTEN__Protein_expression.txt"
data_table <- read.table(file_name, sep = "\t", header = TRUE, quote = "", fill = TRUE)
rppa_table <- read.table(rppa_name, sep = "\t", header = TRUE)

# remove NA
rppa_table <- rppa_table[!is.na(data_table$Overall.Survival..Months.),]
data_table <- data_table[!is.na(data_table$Overall.Survival..Months.),]

# set max time
max_time <- 100
rppa_table <- rppa_table[data_table$Overall.Survival..Months. <= max_time,]
data_table <- data_table[data_table$Overall.Survival..Months. <= max_time,]

# split data into low and high expression
expres <- rppa_table$PTEN..Protein.expression..RPPA.
cut <- quantile(expres, c(0.5, 0.5))
low_expres <- data_table[rppa_table$PTEN..Protein.expression..RPPA. <= cut[1],]
high_expres <- data_table[rppa_table$PTEN..Protein.expression..RPPA. >= cut[2],]

# status: 1=censored 2=dead
low_status <- rep(1, NROW(low_expres))
high_status <- rep(1, NROW(high_expres))
low_status[low_expres$Overall.Survival.Status == "1:DECEASED"] <- 2
high_status[high_expres$Overall.Survival.Status == "1:DECEASED"] <- 2
low_data = list("time"=low_expres$Overall.Survival..Months., "status"=low_status)
high_data = list("time"=high_expres$Overall.Survival..Months., "status"=high_status)

total_data = list("time"=c(low_expres$Overall.Survival..Months.,high_expres$Overall.Survival..Months.),
                  "status"=c(low_status,high_status),
                  "expres"=c(rep(1, NROW(low_status)), rep(2, NROW(high_status))))

# survival analysis
plot(survfit(Surv(time, status)~1, data=low_data),
           xlab = "Months",
           ylab = "Overall survival probability",
           col = "blue")
lines(survfit(Surv(time, status)~1, data=high_data),
           xlab = "Months",
           ylab = "Overall survival probability",
           col = "red")
plot(survfit(Surv(time, status)~expres, data=total_data),
     xlab = "Months",
     ylab = "Overall survival probability")
surv_pvalue(survfit(Surv(time, status)~expres, data=total_data))

