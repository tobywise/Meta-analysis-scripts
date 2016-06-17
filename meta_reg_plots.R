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
  print(sdm_table)
  
  col_names <- colnames(sdm_table)
  
  for (i in 1:length(files)) {
    if(grepl(paste("extract_", analysis_name, ".+.txt",sep = ""), files[i])) {  # only use the files containing extracted values
      name <- gsub(".txt", "", files[i])  # remove .txt part for naming
      print(name)
      data <- read.table(paste(dir, files[i], sep = ""), header = TRUE)[1:3]  # read first three columns of the table
      colnames(data)[1] <- 'study'
      colnames(data)[2] <- paste(name, "_estimate", sep="")
      colnames(data)[3]<- paste(name, "_variance", sep="")
      sdm_table <- merge(sdm_table, data, by='study')
    }
  }
  
  data <- subset(sdm_table, grepl(regex_selector, data$map, perl = TRUE))  # need to select only combined studies here
  
  return(sdm_table)
}


reg_plot <- function(data, extracted_metareg, predictor, save_plots=TRUE, bar=FALSE, invert=FALSE, return_meta=FALSE){
  # Function to create meta-regression plots with point sizes from weights
  # and regression line from meta-regression
  #
  # Arguments:
  # data = data
  # extracted_metareg = name of the meta-regression as a string
  # predictor = predictor variable for meta-regression
  # save_plots = option to save plots as pdfs
  # bar = gives bar plot instead of scatter
  # invert = inverts the y axis so positive effect sizes become negative
  # return_meta = returns a meta analysis on the extracted values rather than the plot
  # 
  # Returns:
  # Outputs plot and saves as a pdf
  #
  
  require(ggplot2)
  require(metafor)
  
  estimates <- colnames(data)[grepl(paste(extracted_metareg, ".+estimate", sep=''), colnames(data))] 
  variances <- colnames(data)[grepl(paste(extracted_metareg, ".+variance", sep=''), colnames(data))] 
  
  data <- subset(data, data[,predictor]!='NA')  # exclude any studies where the predictor value is NA
  print(data[,predictor])
  
  out_list <- list()
  meta_out_list <- list()
  
  for (k in 1:length(estimates)) {
    estimate = estimates[k]
    variance = variances[k]
    
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
      plot_name <- paste(estimate, predictor, sep ="_")
      ggobj <- ggplot(data=data, aes_string(x = predictor, y = estimate)) + geom_point(size = (weights*0.5+1), colour="Gray") + 
        geom_point(shape = 1, size = (weights*0.5+1), colour="black") + theme_bw() +
        geom_line(data = meta_data, aes(x = number, y = preds.pred), size = 0.5) + xlab(predictor)  + ylab('Estimate') +
        scale_x_continuous(limits = c(pred_min, pred_max)) + scale_y_continuous(limits = c(est_min, est_max)) + ggtitle(plot_name)
      if (save_plots == TRUE) {
        ggsave(paste(getwd(), '/', plot_name, '.pdf', sep=''))
      }
      out_list[[plot_name]]<- ggobj
      meta_out_list[[plot_name]] <- meta
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
      plot_name <- paste(estimate, predictor, sep ="_")
      ggobj <- ggplot() + 
        geom_jitter(aes_string(x=predictor, y=estimate), data = data, colour = "Gray", position = position_jitter(width = 0.05), 
                    size = c(weights(meta1)*0.5+1, weights(meta2)*0.5+1)) +
        geom_crossbar(data=meta_b, aes(x = Predictor_var, y = Est, ymin = Est, ymax = Est)) + 
        xlab(predictor)  + ylab(estimate) + theme_bw() + ggtitle(plot_name)
      if (save_plots == TRUE) {
        ggsave(paste(getwd(), '/', plot_name, '.pdf', sep=''))
      }
      out_list[[plot_name]]<- ggobj
      meta_out_list[[plot_name]] <- meta_b
    }
  }
  #ggsave(file=paste(predictor, estimate, ".pdf", sep = "_"), width = 8, height = 6)  # save as pdf
  if (return_meta==TRUE){
    return(meta_out_list)
  }
  else {
    return(out_list)
  }
  #return(meta_b)
}

# Example usage

# Get all extracted values in a given folder and add them to the SDM table
# This adds columns to the SDM table for each meta-regression that represent the effect size
# and variance for each study

sdm_table_metaregressions <- add_extracted("C:/Users/Toby/Documents/bd_metaregs/")

# Create plots from the new SDM table using the lithium meta-regression
# The first argument is the name of the new SDM table variable, the second is the name of 
# the meta-regression as it was called in SDM, the third is the name of the predictor
# variable in the SDM table.
#
# This will create a list of ggplot objects that can be referenced using the $ symbol
# e.g in this example, lithium_plots$plot_name_1 would give you the first plot
# You can check the names of the plot by just calling e.g. lithium_plots
#
# This function also saves the plots as pdfs in the current working directory
# (you can check which directory this is using getwd(), and can change it using
# setwd('path/to/directory'))

lithium_plots <- reg_plot(sdm_table_metaregressions, 'lithium', 'LithiumPercent')

# This function has a couple of extra options, including the ability to do bar plots for
# categorical predictors if needed.



