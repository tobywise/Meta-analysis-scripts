require(ggplot2)
require(metafor)

# extract values and from meta-regressions and adds them to sdm table

dir <- "C:/Users/k1327409/Dropbox/PhD/MDD BD Meta-analysis/Results/BD Results/"
files <- list.files(path=dir)  # gets files in directory

sdm_table <- read.table(paste(dir,'sdm_table.txt', sep = ""), header = TRUE)  # load sdm table
sdm_table <- subset(sdm_table, grepl("et_al_.+[^ab]_sMRI", sdm_table$study))

col_names <- colnames(sdm_table)

for (i in 1:length(files)) {
  if(grepl("extract.+.txt", files[i])) {  # only use the files containing extracted values
    name <- gsub(".txt", "", files[i])  # remove .txt part for naming
    print(name)
    data <- read.table(paste(dir, files[i], sep = ""), header = TRUE)  # read the table
    data <- subset(data, grepl("et_al_.+[^ab]_sMRI", data$map))  # need to select only combined studies here
    est_name <- paste(name, "_estimate", sep="")
    var_name <- paste(name, "_variance", sep="")
    col_names <- c(col_names, est_name, var_name)  # add cluster estimate & variance names to column names
    sdm_table <- cbind(sdm_table, data$estimate)  # append estimate and variance to data
    sdm_table <- cbind(sdm_table, data$variance)
  }
}

colnames(sdm_table) <- col_names  # name the columns properly


# function to create meta-regression plots
reg_plot <- function(data, estimate, variance, predictor){
  # Function to create meta-regression plots with point sizes from weights
  # and regression line from meta-regression
  #
  # Arguments:
  # data = data
  # estimate = estimate extracted from coordinate
  # variance = variance extracted from coordinate
  # predictor = predictor variable for meta-regression
  # 
  # Returns:
  # Outputs plot and saves as a pdf
  #
  data <- subset(data, data[,predictor]!='NA')  # exclude any studies where the predictor value is NA
  meta <- rma(data[,estimate], data[,variance], mods=data[,predictor], data = data)  # do meta-regression
  preds <- predict(meta, newmods = c(0:100))  # predict for 0:100 using meta-regression
  number <- c(0:100)  # create x-axis
  meta_data <- data.frame(number, preds$pred)  # create data frame for predictions
  weights <- weights(meta)  # get weights from meta-regression, use for size in plot
  ggobj <- ggplot(data=data, aes_string(x = predictor, y = estimate)) + geom_point(size = (weights*1.2)) + theme_bw() +
    geom_line(data = meta_data, aes(x = number, y = preds.pred), size = 0.5) + xlab(predictor)  + ylab(estimate)
  ggsave(file=paste(predictor, estimate, ".pdf", sep = "_"), width = 8, height = 6)  # save as pdf
  return(ggobj)
}



# Run this function for all meta-regression outputs
ests <- colnames(sdm_table)[grepl("extract.+estimate", colnames(sdm_table))]
vars <- colnames(sdm_table)[grepl("extract.+variance", colnames(sdm_table))]

for (i in 1:length(ests)) {  # assumes ests and vars are in the same order
  print(ests[i])
  print(vars[i])
  assign(paste(ests[i], "_meds_plot", sep =""), reg_plot(sdm_table, ests[i], vars[i], 'MedicationPercent'))
}


