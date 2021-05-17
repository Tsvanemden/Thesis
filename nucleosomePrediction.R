
# install the package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NuPoP")

# load the package
library("NuPoP")


# perform the prediction on, in this case Fastafile.txt
results <- predNuPoP("Fastafile.txt",species=9,model=4)


#make the plots
plotNuPoP(readNuPoP('Fastafile.txt_Prediction4.txt',1,20000))

