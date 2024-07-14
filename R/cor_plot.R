#' @title Correlation Plot with Significance Points
#' @name cor_plot
#' @description This function generates a correlation plot between two datasets, displaying correlation coefficients as a heatmap and significant correlations as scatter points.
#'
#' @param data1 A data frame or matrix representing the first dataset.
#' @param data2 A data frame or matrix representing the second dataset.
#' @param sig.name1 A character string specifying the name of the first dataset (default: "value1").
#' @param sig.name2 A character string specifying the name of the second dataset (default: "value2").
#' @param cutoff.pvalue The significance threshold for correlation (default: 0.05).
#' @param color A vector of two colors for the heatmap gradient (default: c("#62CCC9", "#FF9999")).
#'
#' @return A ggplot object displaying the correlation heatmap and scatter points.
#'
#' @details The function computes correlation coefficients between corresponding columns in the two datasets and identifies significant correlations based on p-values.
#'
#' @examples
#' # Read all data into memory
#' data(data.path, package = "iPRISM")
#' data(data.cell, package = "iPRISM")
#' # Draw the plot
#' cor_plot(data1 = data.path,data2 = data.cell,sig.name1 = "path",sig.name2 = "cell")
#'
#' @importFrom ggplot2 ggplot aes_string ylab xlab theme
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom Hmisc rcorr
#' @importFrom tidyr pivot_longer
#' @keywords correlation plot heatmap scatter ggplot
#'
#' @export
cor_plot <- function(data1,
                     data2,
                     sig.name1 = "value1",
                     sig.name2 = "value2",
                     cutoff.pvalue = 0.05,
                     color = c("#62CCC9","#FF9999")){
  data_both <- cbind(data1,data2[rownames(data1),])
  rp_data <- rcorr(data_both)
  data_r <- as.data.frame(rp_data$r[-(1:ncol(data1)),-((ncol(data1) + 1):ncol(data_both))])
  data_p <- rp_data$P
  data_r <- as.data.frame(cbind(rownames(data_r),data_r))
  colnames(data_r)[1] <- "value1"
  rp_data <- pivot_longer(data_r,cols = colnames(data_r)[-1],names_to = "value2",values_to = "Rvalue")
  rp_data$Pvalue <- 0
  rp_data <- as.data.frame(rp_data)
  for(i in 1:nrow(rp_data)){
    rp_data[i,4] <- data_p[rp_data[i,1],rp_data[i,2]]
  }
  rp_data$shape <- ifelse(rp_data$Pvalue < cutoff.pvalue,"A","B")
  rp_data$log10p <- ifelse(rp_data$Pvalue < cutoff.pvalue,-log10(rp_data$Pvalue),10)
  rp_data$log10p <- as.numeric(rp_data$log10p)
  ggplot(data = rp_data) +
    geom_tile(data = rp_data,aes_string(x = "value1",y = "value2",fill = "Rvalue"),size = 2,color = "white") +
    geom_point(aes_string(x = "value1",y = "value2",shape = "shape",size = "log10p"),color = "black") +
    scale_shape_manual(values = c(19,4)) +
    theme_bw() +
    scale_fill_gradient2(low = color[1],high = color[2],midpoint = 0) +
    ylab(sig.name2)+
    xlab(sig.name1)+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
}
