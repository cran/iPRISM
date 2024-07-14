## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Load the package
library(iPRISM)

## -----------------------------------------------------------------------------
# Read cell line and pathway information
data(data.path, package = "iPRISM")
data(data.cell, package = "iPRISM")

# Draw the plot
cor_plot(data1 = data.path, data2 = data.cell, sig.name1 = "path", sig.name2 = "cell")

## -----------------------------------------------------------------------------
# Load example data
data(Seeds, package = "iPRISM")
data(ppi, package = "iPRISM")
data(path_list, package = "iPRISM")

# Shrink pathway list to the top 2 pathways
path_list <- path_list[1:5]

# Perform GSEA
result <- get_gsea_path(seed = Seeds, network = ppi, pathlist = path_list, gsea.nperm = 100)
print(result)

## -----------------------------------------------------------------------------
# Load example data
data(data_sig, package = "iPRISM")

# Fit logistic regression model
b <- get_logiModel(data.sig = data_sig, pred.value = pred_value, step = TRUE)
summary(b)

