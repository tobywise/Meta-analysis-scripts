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
  sdm_table <- subset(sdm_table, grepl(regex_selector, sdm_table$study, perl = TRUE))
  #sdm_table$PercentMale <- sdm_table$PercentMale*100  # what is this doing here?!
  print(sdm_table)
  
  col_names <- colnames(sdm_table)
  
  for (i in 1:length(files)) {
    if(grepl(paste("extract_", analysis_name, ".+.txt",sep = ""), files[i])) {  # only use the files containing extracted values
      name <- gsub(".txt", "", files[i])  # remove .txt part for naming
      print(name)
      data <- read.table(paste(dir, files[i], sep = ""), header = TRUE)  # read the table
      data <- head(data, -3)  # remove last 3 rows, which give overall results
      data <- subset(data, grepl(regex_selector, data$map, perl = TRUE))  # need to select only combined studies here
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


reg_plot <- function(data, estimate, variance, predictor,bar=FALSE, invert=FALSE){
  # Function to create meta-regression plots with point sizes from weights
  # and regression line from meta-regression
  #
  # Arguments:
  # data = data
  # estimate = estimate extracted from coordinate
  # variance = variance extracted from coordinate
  # predictor = predictor variable for meta-regression
  # bar = gives bar plot instead of scatter
  # invert = inverts the y axis so positive effect sizes become negative
  # 
  # Returns:
  # Outputs plot and saves as a pdf
  #
  data <- subset(data, data[,predictor]!='NA')  # exclude any studies where the predictor value is NA
  if(bar==FALSE) {
    meta <- rma(data[,estimate], data[,variance], mods=data[,predictor], data = data)  # do meta-regression
    preds <- predict(meta, newmods = c(0:100))  # predict for 0:100 using meta-regression
    number <- c(0:100)  # create x-axis
    print(data[,estimate])
    if(invert==TRUE){
      data[,estimate] <- -data[,estimate]
      preds$pred <- -preds$pred
    }
    print(data[,estimate])
    meta_data <- data.frame(number, preds$pred)  # create data frame for predictions
    weights <- weights(meta)  # get weights from meta-regression, use for size in plot    
    pred_max = max(data[,predictor]) + (max(data[,predictor]))/10
    pred_min = min(data[,predictor]) - (max(data[,predictor]))/10
    est_max = max(data[,estimate]) + (max(data[,estimate]))/10
    est_min = min(data[,estimate]) - (max(data[,estimate]))/10
    ggobj <- ggplot(data=data, aes_string(x = predictor, y = estimate)) + geom_point(size = (weights*0.5+1)) + theme_bw() +
      geom_line(data = meta_data, aes(x = number, y = preds.pred), size = 0.5) + xlab(predictor)  + ylab(estimate) +
      scale_x_continuous(limits = c(0, pred_max)) + scale_y_continuous(limits = c(est_min, est_max))
  }
  
  else {
    if(invert==TRUE){
      data[,estimate] <- -data[,estimate]
    }
    meta1 <- rma(data[,estimate][data[, predictor]==0], data[,variance][data[, predictor]==0], data = data)
    meta2 <- rma(data[,estimate][data[, predictor]==1], data[,variance][data[, predictor]==1], data = data)
    meta_b <- data.frame(c(meta1[['b']], meta2[['b']]), c("0","1"), c(meta1[['ci.lb']], meta2[['ci.lb']]), c(meta1[['ci.ub']], meta2[['ci.ub']]),
                         c(meta1[['se']], meta2[['se']]))
                                                                      
    colnames(meta_b) <- c("Est", "Predictor_var", "LB", "UB", "SE")
    print(meta_b)
    group_stats <- matrix(, nrow = length(unique(data[,predictor])), ncol = 3)
    group_stats <- data.frame(group_stats)
    count = 1
    for (i in unique(data[,predictor])){
      vals <- data[,estimate][data[,predictor]==i]
      group_mean <- mean(vals)
      group_sd <- sd(vals)
      stats_vec <- c(i, group_mean, group_sd)
      group_stats[count,] <- stats_vec
      group_stats[,2] <- as.numeric(group_stats[,2]) #awkward solution to discrete problem
      group_stats[,3] <- as.numeric(group_stats[,3])
      count <- count +1
    }
    data[,predictor] <- as.factor(data[,predictor])
    meta_b$Predictor_var <- as.factor(meta_b$Predictor_var)
    colnames(group_stats) <- c("Predictor_var", "Mean", "SD")
    ggobj <- ggplot() + 
      geom_jitter(aes_string(x=predictor, y=estimate), data = data, colour = I("dodgerblue4"), position = position_jitter(width = 0.05), 
                  size = c(weights(meta1)*0.5+1, weights(meta2)*0.5+1)) +
      geom_crossbar(data=meta_b, aes(x = Predictor_var, y = Est, ymin = Est, ymax = Est)) + 
      xlab(predictor)  + ylab(estimate) + theme_bw()
  }
  
  #ggsave(file=paste(predictor, estimate, ".pdf", sep = "_"), width = 8, height = 6)  # save as pdf
  #return(ggobj)
  return(meta_b)
}

#Example usage

# Get extracted values for all coordinates in a folder
sdm_table <- add_extracted("C:/Users/k1327409/Documents/VBShare/24_03/", "", "^(combined_t_)?[A-Za-z]+_?[A-Za-z]+_et_al_.+(?<![012456789][ab])_sMRI")

# Run plotting function for all meta-regression outputs
ests <- colnames(sdm_table)[grepl("extract_2303_BD_dep_euth.+estimate", colnames(sdm_table))]  # get estimates
vars <- colnames(sdm_table)[grepl("extract_2303_BD_dep_euth.+variance", colnames(sdm_table))]  # get variances


for (i in 1:length(ests)) {  # assumes ests and vars are in the same order
  print(ests[i])
  print(vars[i])
  assign(paste(ests[i], "_mood_state_plot", sep =""), reg_plot(sdm_table_dep_euth, ests[i], vars[i], 'State'))
}


antipsychotics_multi <- multiplot(extract_2303_BD_antipsychotic_metareg_coords_1_estimate_antipsychotic_plot + 
                                    ylab("C1 Est"), extract_2303_BD_antipsychotic_metareg_coords_2_estimate_antipsychotic_plot + 
                                    ylab("C2 Est"), extract_2303_BD_antipsychotic_metareg_coords_3_estimate_antipsychotic_plot + 
                                    ylab("C3 Est"),extract_2303_BD_antipsychotic_metareg_coords_4_estimate_antipsychotic_plot + 
                                    ylab("C4 Est"),extract_2303_BD_antipsychotic_metareg_coords_5_estimate_antipsychotic_plot + 
                                    ylab("C5 Est"), cols=2)
lithium_multi <- multiplot(extract_2303_BD_lithium_metareg_coords_1_estimate_lithium_plot + ylab("C1 Est"), extract_2303_BD_lithium_metareg_coords_2_estimate_lithium_plot + 
                                    ylab("C2 Est"), cols=2)
sex_multi <- multiplot(extract_2303_BD_sex_metareg_coords_1_estimate_sex_plot + 
                                    ylab("C1 Est"), extract_2303_BD_sex_metareg_coords_2_estimate_sex_plot + 
                                    ylab("C2 Est"), extract_2303_BD_sex_metareg_coords_3_estimate_sex_plot + 
                                    ylab("C3 Est"),extract_2303_BD_sex_metareg_coords_4_estimate_sex_plot + 
                                    ylab("C4 Est"), cols=2)
state_multi <- multiplot(extract_2303_BD_dep_euth_metareg_coords_1_estimate_mood_state_plot + 
                         ylab("C1 Est"), extract_2303_BD_dep_euth_metareg_coords_2_estimate_mood_state_plot + 
                         ylab("C2 Est"), extract_2303_BD_dep_euth_metareg_coords_3_estimate_mood_state_plot + 
                         ylab("C3 Est"),extract_2303_BD_dep_euth_metareg_coords_4_estimate_mood_state_plot + 
                         ylab("C4 Est"), cols=2)


# MDD
# Get extracted values for all coordinates in a folder
sdm_table_mdd <- add_extracted("C:/Users/k1327409/Documents/VBShare/MDD_sMRI/Meta_regressions/Extracted/", "", "^(combined_t_)?[A-Za-z]+[_-]?[A-Za-z]+_et_al_.+[^abc]_sMRI")

# Run plotting function for all meta-regression outputs
ests_mdd <- colnames(sdm_table_mdd)[grepl("extract_10_03_MDD_Sex.+estimate", colnames(sdm_table_mdd))]  # get estimates
vars_mdd <- colnames(sdm_table_mdd)[grepl("extract_10_03_MDD_Sex.+variance", colnames(sdm_table_mdd))]  # get variances

for (i in 1:length(ests_mdd)) {  # assumes ests and vars are in the same order
  print(ests_mdd[i])
  print(vars_mdd[i])
  print(paste(ests_mdd[i], "_MDD_Sex_plot", sep =""))
  assign(paste(ests_mdd[i], "_MDD_Sex_plot", sep =""), reg_plot(sdm_table_mdd, ests_mdd[i], vars_mdd[i], 'PercentMale', invert=TRUE))
}

mdd_severity_multi <- multiplot(extract_10_03_MDD_Severity_coords_1_estimate_MDD_Severity_plot + 
                                    ylab("C1 Est"), extract_10_03_MDD_Severity_coords_2_estimate_MDD_Severity_plot + 
                                    ylab("C2 Est"), extract_10_03_MDD_Severity_coords_3_estimate_MDD_Severity_plot + 
                                    ylab("C3 Est"),extract_10_03_MDD_Severity_coords_4_estimate_MDD_Severity_plot, cols=2)
mdd_antidepressants_multi <- multiplot(extract_10_03_MDD_Antidepressants_coords_1_estimate_MDD_Antidepressants_plot + ylab("C1 Est"), 
                                       extract_10_03_MDD_Antidepressants_coords_2_estimate_MDD_Antidepressants_plot + ylab("C2 Est"), cols=2)
mdd_sex_multi <- multiplot(extract_MDD_sex_coords_1_estimate_MDD_Sex_plot + 
                         ylab("C1 Est"), extract_MDD_sex_coords_2_estimate_MDD_Sex_plot + 
                         ylab("C2 Est"), extract_MDD_sex_coords_3_estimate_MDD_Sex_plot + 
                         ylab("C3 Est"), cols=2)

extract_10_03_MDD_Severity_coords_2_estimate_MDD_Severity_plot + 
  ylab("SDM Z") + xlab("HAM-D Score")

extract_10_03_MDD_Sex_coords_2_estimate_MDD_Sex_plot + 
  ylab("SDM Z") + xlab("Percentage of Male Patients")