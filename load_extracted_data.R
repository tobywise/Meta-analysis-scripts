# load extracted info for clusters and do funnel plots & egger's test

require(metafor)

dir <- "C:/Users/k1327409/Dropbox/PhD/MDD BD Meta-analysis/Results/BD Results/"
files <- list.files(path=dir) 

for (i in 1:length(files)) {
  if(grepl("extract.+.txt", files[i])) {
    name <- gsub(".txt", "", files[i])
    print(name)
    data <- read.table(paste(dir, files[i], sep = ""), header = TRUE)
    data <- subset(data, grepl("et_al_.+[^ab]_sMRI", data$map))
    assign(name, data)
    ma <- rma(data$estimate, data$variance)
    funnel(ma) 
    assign(paste(name, "_funnel", sep= ""), recordPlot())
    assign(paste(name, "_ma", sep= ""), ma)
    assign(paste(name, "_egger", sep= ""), regtest(ma))
    print(regtest(ma))
  }
}