library("pegas")
library("vcfR")

setwd('../01_simulating-test-data/data/')

write.table('vcf_file n_missing pi theta tajimad n_sites', 
            file = 'pegas_summary.txt', col.names = FALSE, row.names = FALSE,
            quote = FALSE)

VCF_files = list.files(getwd(), recursive = TRUE)

for (VCF_file in VCF_files) {

VCF = read.vcfR(VCF_file)
VCF@fix[,4] <- rep("A", length(VCF@fix[,4]))
VCF@fix[,5] <- rep("T", length(VCF@fix[,5]))
DNABIN = vcfR2DNAbin(VCF)

Pi = nuc.div(DNABIN)
Theta = theta.s(DNABIN)
TajimaD = tajima.test(DNABIN)[1]

if(grepl("missing_genos", VCF_file)){
  
  n_missing <- VCF_file %>% gsub(".*invar_missing_genos=", "", .) %>% 
    gsub(".vcf.gz", "", .) %>%
    as.numeric
  
  n_missing <- 10000*n_missing
  
  
} else{
  
  n_missing <- VCF_file %>% gsub(".*invar.missing_", "", .) %>% 
    gsub(".vcf", "", .) %>%
    as.numeric
  
  n_missing <- 10000 - n_missing
  
}

write.table(paste(VCF_file, n_missing, Pi, Theta, TajimaD, 10000), 
            file = 'pegas_summary.txt', col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = TRUE)
}

