#' @title Fit Logistic Regression Model
#' @name get_logiModel
#' @description This function fits a logistic regression model to the given data.
#'
#' @param data.sig A data frame where each row is a sample and each column is a pathway.
#' @param pred.value A numeric vector representing the response variable.
#' @param levels A character vector specifying the levels of the response variable (default: c("R", "N")).
#' @param step Logical. If TRUE, perform stepwise model selection (default: TRUE).
#'
#' @return A fitted logistic regression model.
#'
#' @details The function converts the response variable to a factor with specified levels and fits a logistic regression model using the glm function.
#'
#' @examples
#' data(data_sig, package = "iPRISM")
#' \donttest{
#' b <- get_logiModel(data.sig = data_sig, pred.value = pred_value, step = TRUE)
#' summary(b)
#' }
#' @importFrom stats glm
#' @keywords logistic regression model classification
#'
#' @export
get_logiModel <- function(data.sig,
                          pred.value,
                          levels = c("R","N"),
                          step = TRUE){
  data.sig <- as.data.frame(data.sig)
  data.sig$pred <- pred.value[rownames(data.sig)]
  data.sig$pred <- factor(data.sig$pred,levels = levels,ordered = T)
  model <- glm(pred ~ .,data = data.sig,family = "binomial")
  if(step == TRUE){
    model <- step(object = model,trace = 0)
  }
  return(model)
}
