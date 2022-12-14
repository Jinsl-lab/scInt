# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

norm_S_batches <- function(x, p, i, ncol, groups, N) {
    .Call('_scInt_norm_S_batches', PACKAGE = 'scInt', x, p, i, ncol, groups, N)
}

upper_tri <- function(x, p, i, ncol) {
    .Call('_scInt_upper_tri', PACKAGE = 'scInt', x, p, i, ncol)
}

filter_T <- function(x, p, i, ncol, T_th, A, B) {
    .Call('_scInt_filter_T', PACKAGE = 'scInt', x, p, i, ncol, T_th, A, B)
}

eigenMapMatMult <- function(A, B) {
    .Call('_scInt_eigenMapMatMult', PACKAGE = 'scInt', A, B)
}

eigenMapMatcrossprod <- function(A, B) {
    .Call('_scInt_eigenMapMatcrossprod', PACKAGE = 'scInt', A, B)
}

eigenMapMattcrossprod <- function(A, B) {
    .Call('_scInt_eigenMapMattcrossprod', PACKAGE = 'scInt', A, B)
}

