tumours <- list.files(path = "TCGA_all_data")
tumours <- tumours[tumours != "COADREAD"]
tumours <- tumours[tumours != "STES"]
tumours <- tumours[tumours != "KIPAN"]
tumours <- tumours[tumours != "GBMLGG"]

patients_with_CNV_info <- list()
for(tumour in tumours){
  TCGA_files <- read.table(paste("TCGA_all_data/", tumour, "/gdac.broadinstitute.org_", tumour, ".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/", tumour, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep=""),
                           header=TRUE)  #File from TCGA

  #Select tumour samples
  samples <- TCGA_files$Sample
  sample_codes <- substring(samples, 14, 16)
  
  #Sample types available: "10A" "01A" "11A" "11B" "01B" "10B" "06A"
  
  tumour_samples <- which(sample_codes %in% c("01A", "01B", "01C", "01D"))
  
  TCGA_files <- TCGA_files[c(tumour_samples),]
  patients <- unique(substr(TCGA_files$Sample, 9, 12))
  patients_with_CNV_info[[tumour]] <- patients
  print(tumour)
}
save(patients_with_CNV_info, file="patients_with_CNV_info.Rdata")

sapply(patients_with_CNV_info, length)
