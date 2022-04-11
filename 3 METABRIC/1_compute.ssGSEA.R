#' Purpose:
#' To compute ssGSEAs (ie. Find Enrichment of Signatures in BC expression profiles)

# Call required packages into environment
library(dplyr)
library(GSVA)
library(readxl)
library(foreach)
library(ggfortify)

# Get today's date in a variable
todays_date <- format(Sys.Date(), "%Y%m%d")

#loading the signatures Function
loadSignatures = function(file){
  sign_file = read_xlsx(file)
  sign_list = lapply(sign_file, function(x) x[!is.na(x)])
  return(sign_list)
}

# Load RNA seq METABRIC data
load(file="metabricExpression.RData") # Note: too large to upload
metabricExpression = metabricExpression %>% select(Symbol, everything())

# Function to filter and clean metabricExpression.RData object
CleanRnaSeqFunc <- function(rnaSeq) {
  rnaSeq = metabricExpression[!duplicated(metabricExpression$Symbol),]
  rownames(rnaSeq) <- rnaSeq$Symbol
  rnaSeq = rnaSeq %>% select(
    -CLID,
    -Symbol
  )
  # Remove duplicated rows
  rnaSeq <- rnaSeq[!duplicated(colnames(rnaSeq)),]
  rnaSeq.clean = as.matrix(as.data.frame(rnaSeq))
}
rnaSeq.clean = CleanRnaSeqFunc(metabricExpression)

# Get METABRIC clinical annotations
load("metabricClinical.RData")

# Get file name
signature_file <- list.files(pattern="xlsx")

# Create output folder 
output_folder <- sprintf("%s Analysis", gsub("\\.xlsx", "" , signature_file))
dir.create(output_folder)

# Import signatures
signatures = loadSignatures(signature_file)  

# Compute ssGSEA matrix
ssGSEA = gsva(
  rnaSeq.clean, 
  signatures, 
  method="ssgsea", 
  kcdf="Gaussian")

# Add the individual patient id as a separate column
ssGSEA.ind <- t(ssGSEA)
rownames(ssGSEA.ind) <- gsub("\\.", "-", rownames(ssGSEA.ind))
ssGSEA.ind.df <- as.data.frame(ssGSEA.ind) %>% mutate(
  Patient.ID = as.character(rownames(ssGSEA.ind))
)

# Merge exp data with clinical annotatons
ssGSEA.subtypes.file <- merge(
  ssGSEA.ind.df,
  metabricClinical,
  by="Patient.ID",
  x.all = T
)
rownames(ssGSEA.subtypes.file) <- ssGSEA.subtypes.file$Patient.ID

# Filter out NA subtypes
ssGSEAs.subt <- ssGSEA.subtypes.file %>% filter(
  !Pam50...Claudin.low.subtype == "NA"
) %>% select(
  Pam50...Claudin.low.subtype,
  everything(),
  -Sample.ID
) 

# Rename new ssGSEA matrix rows
rownames(ssGSEAs.subt) <- ssGSEAs.subt$Patient.ID
ssGSEAs.subt.clinical <- ssGSEAs.subt %>% select(
  -Patient.ID
)
ssGSEAs <- ssGSEAs.subt.clinical[1:5]

# Save ssGSEAs to files
#ssgsea_csv <-  "ssGSEA_3sign%s.csv"
# write.csv(ssGSEAs, file = ssgsea_csv , row.names=T)
# BC subtypes of interest
subtypes <- c("Basal", "claudin-low", "LumA", "LumB", "Her2")

# Create ssGSEAs file by subtypes
ssGSEA.subtype.list  = foreach (subtype = 1:length(subtypes)) %do% {
  ssGSEAs.subt.func <- ssGSEAs %>% mutate(
    ind = rownames(ssGSEAs)
  )
  subtype.file =
    ssGSEAs.subt.func %>% filter(
      Pam50...Claudin.low.subtype == subtypes[subtype]
    )
  rownames(subtype.file) <- subtype.file$ind
  subtype.file = subtype.file %>% select(
    -ind
  )
  # ssGSEAs.subtype.filename = sprintf(
  #   "%s/ssGSEA_3sign_%s_%s.csv", metabric.subtype.Path[subtype], "metabric",
  #   subtypes[subtype])
  # write.csv(subtype.file, ssGSEAs.subtype.filename)
  return(subtype.file)
}
# Rename list elements
names(ssGSEA.subtype.list) <- subtypes

# Save R data object as rds file
subtype_ssgsea_rds <-  sprintf("%s/ssGSEA.subtype.list.3sign.RDS", output_folder)
saveRDS(ssGSEA.subtype.list,file =subtype_ssgsea_rds)

