## Vectorize then parallelize with respect to an argument 
## See ?parallel::pvec
pvec2 = function(v, FUN, arg=NULL, mc.cores = NULL, mc.silent = FALSE, mc.cleanup = TRUE, ...){
  if(is.null(mc.cores)){
    mc.cores = parallel::detectCores() - 1
  }
  arg = ifelse(is.null(arg), formalArgs(args(FUN))[1], arg)
  parallel::pvec(v, Vectorize(FUN, arg), ..., mc.cores = mc.cores, mc.silent = mc.silent, mc.cleanup = mc.cleanup)
}

## Get the next possible floating point number from x in the direction of y
Rcpp::cppFunction(
  code = 
'double nextAfter(double x, double y){
  return std::nextafter(x, y);
}')

## Create a latex table with a caption using gt (because gt
## doesn't gracefully support this)
## Adapted from https://github.com/rstudio/gt/issues/818#issuecomment-1287060966
as_latex_with_caption <- function(gt, chunk_label, afterpage = FALSE, spacing=NULL) {
  if(!knitr::is_latex_output()) return(gt)
  cap_text = gt$`_options`$value[which(gt$`_options`$parameter == "table_caption")][[1]]
  gt <- gt::as_latex(gt)
  caption <- glue::glue("\\caption{{{cap_text}\\label{{tab:{chunk_label}}}}}\\\\")
  latex <- strsplit(gt[1], split = "\n")[[1]]
  latex <- c(latex[1], caption, latex[-1])
  latex <- paste(latex, collapse = "\n")
  if(!is.null(spacing))
    latex <- glue::glue('\\begin{{spacing}}{{{spacing}}}{latex}\\end{{spacing}}')
  if(afterpage)
    latex <- glue::glue('\\afterpage{{{latex}}}')
  gt[1] <- latex
  return(gt)
}