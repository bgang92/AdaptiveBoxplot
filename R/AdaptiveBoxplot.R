# Helper function for Benjamini-Hochberg procedure (internal)
bh.func<-function(pv, q)
{
  # This function implements the BH procedure
  # the input
  # pv: the p-values
  # q: the FDR level
  # the output
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)
  pvi<-st.pv/1:m
  hps<-rep(0, m)
  if (max(pvi<=(q/m))==0)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}



# Helper function for Holm procedure (internal)
holm.func<-function(pv, q,kfwer=1)
{
  #This function implements the Holm procedure
  # the input
  # pv: the p-values
  # q: the k-FWER level
  # kfwer: k in k-FWER
  # the output
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule

  m=length(pv)
  st.pv<-sort(pv)
  pvi<-st.pv*(m+kfwer-(1:m))
  hps<-rep(0, m)
  if (pvi[1]>q)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k=min(which(pvi>q))
    pk<-st.pv[k-1]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}











# Helper function to calculate BH boxplot statistics (internal)
bh.box=function(x,alpha){
  q <- quantile(x, probs = c(0.25, 0.75))
  xbar =mean(c(q[1],q[2]))
  sd=(q[2]-q[1])/1.35
  zscore=(x-xbar)/sd
  pv=2*(1-pnorm(abs(zscore)))
  bh.res=bh.func(pv,alpha)
  if(bh.res$th==1){
    whis=max(c(abs(min(x)-q[1]), abs(max(x)-q[2])))
  }else{
    x.th=x[which(pv==bh.res$th)]
    whis=min(c(abs(x.th-q[1]), abs(x.th-q[2])))
  }
  result=list(q1=q[1],q3=q[2],med=median(x),lf=q[1]-whis,uf=q[2]+whis,whis=whis)
  return(result)
}

# Helper function to calculate Holm boxplot statistics (internal)
holm.box=function(x,alpha,kfwer=1){

  q <- quantile(x, probs = c(0.25, 0.75))
  xbar =mean(c(q[1],q[2]))
  sd=(q[2]-q[1])/1.35
  zscore=(x-xbar)/sd
  pv=2*(1-pnorm(abs(zscore)))
  holm.res=holm.func(pv,alpha,kfwer)
  if(holm.res$th==1){
    whis=max(c(abs(min(x)-q[1]), abs(max(x)-q[2])))
  }else{
    x.th=x[which(pv==holm.res$th)]
    whis=min(c(abs(x.th-q[1]), abs(x.th-q[2])))
  }
  result=list(q1=q[1],q3=q[2],med=median(x),lf=q[1]-whis,uf=q[2]+whis,whis=whis)
  return(result)
}



#' False Discovery Rate (FDR) Boxplot
#'
#' Generates a boxplot where whisker lengths are determined by the Benjamini-Hochberg
#' procedure to control the False Discovery Rate (FDR), making the outlier
#' detection rule adaptive to sample size and data characteristics.
#'
#' @details This function is a graphical implementation of the p-value pipeline
#' proposed by Gang, Lin, and Tong (2025). It uses robust estimators for the
#' mean and standard deviation based on quartiles to calculate p-values for each
#' observation, then applies the Benjamini-Hochberg (BH) procedure to determine
#' an adaptive p-value threshold for outlier detection. Outliers are points
#' falling beyond the fences defined by this threshold.
#'
#' @param data A numeric vector for a single boxplot, or a data frame for grouped boxplots.
#' @param alpha The target FDR level. Defaults to 0.01.
#' @param group_col A string specifying the name of the grouping column in `data`.
#' @param value_col A string specifying the name of the value column in `data`.
#' @param ... Additional arguments passed to the base \code{\link[graphics]{boxplot}} function.
#'
#' @return A plot is drawn on the current graphics device.
#' @export
#' @seealso \code{\link{holm.boxplot}}
#' @references Gang, B., Lin, H., & Tong, T. (2025). Unifying Boxplots: A Multiple Testing Perspective.
#' @examples
#' # Single group example
#' set.seed(123)
#' data_single <- c(rnorm(50), 10, 12)
#' bh_boxplot(data_single, alpha = 0.05, main = "FDR Boxplot (Single Group)")
#'
#' # Grouped data example
#' data_grouped <- data.frame(
#'   Category = rep(c("A", "B"), each = 100),
#'   Value = c(rnorm(100), rnorm(100, mean = 2, sd = 1.5))
#' )
#' bh_boxplot(data_grouped, group_col = "Category", value_col = "Value")
bh_boxplot=function (data,alpha=0.01, group_col = NULL, value_col = NULL, ...)
{

  if (!is.null(group_col) && !is.null(value_col) && is.data.frame(data)) {
    values <- data[[value_col]]
    group <- factor(data[[group_col]])
  }
  else if (is.null(group_col) && is.null(value_col)) {
    if (is.numeric(data)) {
      values <- data
    }
    else if (is.data.frame(data) && ncol(data) == 1 && is.numeric(data[,
                                                                       1])) {
      values <- data[, 1]
    }
    else {
      stop("For single-group data, provide a numeric vector or a single-column data frame with numeric values.")
    }
    group <- NULL
  }
  else {
    stop("Incorrect data format or arguments.  Provide a numeric vector/data frame for single-group or a data frame with 'group_col' and 'value_col' for grouped data.")
  }
  if (!is.null(group)) {
    BH_ranges <- tapply(values, group, function(x) {
      res=bh.box(x,alpha)
      return(res$whis/IQR(x))

    })
  }
  else {
    n <- length(values)
    res=bh.box(values,alpha)
    BH_range=res$whis/IQR(values)

  }
  if (!is.null(group)) {
    plot(NULL, xlim = c(0.5, length(levels(group)) + 0.5),
         ylim = range(values, na.rm = TRUE), xaxt = "n", xlab = group_col,
         ylab = value_col, main = "BH-type Boxplot")
    axis(1, at = 1:length(levels(group)), labels = levels(group))
    for (i in 1:length(levels(group))) {
      group_name <- levels(group)[i]
      group_data <- values[group == group_name]
      boxplot(group_data, at = i, range = BH_ranges[group_name],
              add = TRUE, names = "", ...)
    }
  }
  else {
    boxplot(values, range = BH_range, outpch = 19,
            main = "BH-type Boxplot", ...)
  }
}



#' Family-Wise Error Rate (FWER) Boxplot
#'
#' Generates a boxplot where whisker lengths are determined by the Holm procedure
#' to control the Family-Wise Error Rate (FWER), providing a conservative yet
#' principled approach to outlier detection.
#'
#' @details This function is a graphical implementation of the p-value pipeline
#' proposed by Gang, Lin, and Tong (2025). It uses robust estimators for the
#' mean and standard deviation based on quartiles to calculate p-values for each
#' observation, then applies the Holm procedure to determine a p-value threshold
#' that controls the FWER. This method is generally more conservative than the
#' FDR boxplot.
#'
#' @param data A numeric vector for a single boxplot, or a data frame for grouped boxplots.
#' @param alpha The target FWER level. Defaults to 0.05.
#' @param kfwer The "k" in k-FWER control. Defaults to 1 for standard FWER.
#' @param group_col A string specifying the name of the grouping column in `data`.
#' @param value_col A string specifying the name of the value column in `data`.
#' @param ... Additional arguments passed to the base \code{\link[graphics]{boxplot}} function.
#'
#' @return A plot is drawn on the current graphics device.
#' @export
#' @seealso \code{\link{bh_boxplot}}
#' @references Gang, B., Lin, H., & Tong, T. (2025). Unifying Boxplots: A Multiple Testing Perspective.
#' @examples
#' # Single group example
#' set.seed(123)
#' data_single <- c(rnorm(50), 10, 12)
#' holm_boxplot(data_single, alpha = 0.05, main = "FWER Boxplot (Single Group)")
#'
#' # Grouped data example
#' data_grouped <- data.frame(
#'   Category = rep(c("A", "B"), each = 100),
#'   Value = c(rnorm(100), rnorm(100, mean = 2, sd = 1.5))
#' )
#' holm_boxplot(data_grouped, group_col = "Category", value_col = "Value")
holm_boxplot=function (data,alpha=0.05,kfwer=1, group_col = NULL, value_col = NULL, ...)
{

  if (!is.null(group_col) && !is.null(value_col) && is.data.frame(data)) {
    values <- data[[value_col]]
    group <- factor(data[[group_col]])
  }
  else if (is.null(group_col) && is.null(value_col)) {
    if (is.numeric(data)) {
      values <- data
    }
    else if (is.data.frame(data) && ncol(data) == 1 && is.numeric(data[,
                                                                       1])) {
      values <- data[, 1]
    }
    else {
      stop("For single-group data, provide a numeric vector or a single-column data frame with numeric values.")
    }
    group <- NULL
  }
  else {
    stop("Incorrect data format or arguments.  Provide a numeric vector/data frame for single-group or a data frame with 'group_col' and 'value_col' for grouped data.")
  }
  if (!is.null(group)) {
    HOLM_ranges <- tapply(values, group, function(x) {
      res=holm.box(x,alpha)
      return(res$whis/IQR(x))

    })
  }
  else {
    n <- length(values)
    res=holm.box(values,alpha)
    HOLM_range=res$whis/IQR(values)

  }
  if (!is.null(group)) {
    plot(NULL, xlim = c(0.5, length(levels(group)) + 0.5),
         ylim = range(values, na.rm = TRUE), xaxt = "n", xlab = group_col,
         ylab = value_col, main = "Holm-type Boxplot")
    axis(1, at = 1:length(levels(group)), labels = levels(group))
    for (i in 1:length(levels(group))) {
      group_name <- levels(group)[i]
      group_data <- values[group == group_name]
      boxplot(group_data, at = i, range = HOLM_ranges[group_name],
              add = TRUE, names = "", ...)
    }
  }
  else {
    boxplot(values, range = HOLM_range, outpch = 19,
            main = "Holm-type Boxplot", ...)
  }
}

