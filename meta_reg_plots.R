require(ggplot2)
require(metafor)

add_extracted <- function(dir, analysis_name, regex_selector){
  # Function to extract values and from meta-regressions and add them to sdm table
  # 
  # Arguments:
  # dir = directory where the results of the extraction are (the .txt files)
  # analysis_name = name of the analysis you want the extracted values for, e.g. if the files are named
  #   extract_MDD_medication_metareg_C1.txt etc., the analysis name would be "MDD_medication_metareg"
  #   This allows you to select only the extracted values for analyses you're interested in. Set to "" 
  #   to load all extracted results
  # regex_selector = regular expression to select particular studies, set to "" if not needed
  # 
  # Returns:
  # sdm_table = sdm_table with additional columns for estimate and variance at each coordinate
  #

  files <- list.files(path=dir)  # gets files in directory
  
  sdm_table <- read.table(paste(dir,'sdm_table.txt', sep = ""), header = TRUE)  # load sdm table
  sdm_table <- subset(sdm_table, grepl(regex_selector, sdm_table$study))
  sdm_table$PercentMale <- sdm_table$PercentMale*100
  
  col_names <- colnames(sdm_table)
  
  for (i in 1:length(files)) {
    if(grepl(paste("extract_", analysis_name, ".+.txt",sep = ""), files[i])) {  # only use the files containing extracted values
      name <- gsub(".txt", "", files[i])  # remove .txt part for naming
      print(name)
      data <- read.table(paste(dir, files[i], sep = ""), header = TRUE)  # read the table
      data <- head(data, -3)  # remove last 3 rows, which give overall results
      data <- subset(data, grepl(regex_selector, data$map))  # need to select only combined studies here
      est_name <- paste(name, "_estimate", sep="")
      var_name <- paste(name, "_variance", sep="")
      col_names <- c(col_names, est_name, var_name)  # add cluster estimate & variance names to column names
      sdm_table <- cbind(sdm_table, data$estimate)  # append estimate and variance to data
      sdm_table <- cbind(sdm_table, data$variance)
    }
  }
  
  colnames(sdm_table) <- col_names  # name the columns properly
  return(sdm_table)
}


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

#Example usage

# Get extracted values for all coordinates in a folder
sdm_table <- add_extracted("C:/Users/k1327409/Dropbox/PhD/MDD BD Meta-analysis/Results/BD Results/1602/", "", "et_al_.+[^ab]_sMRI")

# Run plotting function for all meta-regression outputs
ests <- colnames(sdm_table)[grepl("extract_BD_sex.+estimate", colnames(sdm_table))]  # get estimates
vars <- colnames(sdm_table)[grepl("extract_BD_sex.+variance", colnames(sdm_table))]  # get variances

for (i in 1:length(ests)) {  # assumes ests and vars are in the same order
  print(ests[i])
  print(vars[i])
  assign(paste(ests[i], "_sex_plot", sep =""), reg_plot(sdm_table, ests[i], vars[i], 'PercentMale'))
}

antipsychotics_multi <- multiplot(extract_BD_antipsychotic_metareg_coords_1_estimate_antipsychotic_plot + 
                                    ylab("C1 Est"), extract_BD_antipsychotic_metareg_coords_2_estimate_antipsychotic_plot + 
                                    ylab("C2 Est"), extract_BD_antipsychotic_metareg_coords_3_estimate_antipsychotic_plot + 
                                    ylab("C3 Est"),extract_BD_antipsychotic_metareg_coords_4_estimate_antipsychotic_plot + 
                                    ylab("C4 Est"),extract_BD_antipsychotic_metareg_coords_5_estimate_antipsychotic_plot + 
                                    ylab("C5 Est"),extract_BD_antipsychotic_metareg_coords_6_estimate_antipsychotic_plot + 
                                    ylab("C6 Est"), cols=2)
lithium_multi <- multiplot(extract_BD_lithium_metareg_coords_1_estimate_lithium_plot + ylab("C1 Est"), extract_BD_lithium_metareg_coords_2_estimate_lithium_plot + ylab("C2 Est"), cols=2)
sex_multi <- multiplot(extract_BD_sex_metareg_coords_1_estimate_sex_plot + 
                                    ylab("C1 Est"), extract_BD_sex_metareg_coords_2_estimate_sex_plot + 
                                    ylab("C2 Est"), extract_BD_sex_metareg_coords_3_estimate_sex_plot + 
                                    ylab("C3 Est"),extract_BD_sex_metareg_coords_4_estimate_sex_plot + 
                                    ylab("C4 Est"), cols=2)