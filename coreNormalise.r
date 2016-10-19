setwd("CELfiles/blood/")
# load the oligo library
library(oligo)

# Read in the CEL files in the directory
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)

# You might need to install and load a package for the specific array you are using (this example is mouse gene 2.0 ST)
# It may try to load it automatically, but may fail.  Install & load the library manually if this happens.
library(pd.huex.1.0.st.v2)
core <- rma(affyRaw, target="core")
#probeset <- rma(affyRaw, target="probeset")
#full <- rma(affyRaw, target="full")
#extended <- rma(affyRaw, target="extended")

# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(core,file="core")
#write.exprs(probeset,file="probeset")
#write.exprs(full, file="full")
#write.exprs(extended, file="extended")
