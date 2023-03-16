#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _scInt_eigenMapMatcrossprod(SEXP, SEXP);
extern SEXP _scInt_eigenMapMatMult(SEXP, SEXP);
extern SEXP _scInt_eigenMapMattcrossprod(SEXP, SEXP);
extern SEXP _scInt_filter_T(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _scInt_norm_S_batches(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _scInt_upper_tri(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_scInt_eigenMapMatcrossprod",  (DL_FUNC) &_scInt_eigenMapMatcrossprod,  2},
    {"_scInt_eigenMapMatMult",       (DL_FUNC) &_scInt_eigenMapMatMult,       2},
    {"_scInt_eigenMapMattcrossprod", (DL_FUNC) &_scInt_eigenMapMattcrossprod, 2},
    {"_scInt_filter_T",              (DL_FUNC) &_scInt_filter_T,              7},
    {"_scInt_norm_S_batches",        (DL_FUNC) &_scInt_norm_S_batches,        6},
    {"_scInt_upper_tri",             (DL_FUNC) &_scInt_upper_tri,             4},
    {NULL, NULL, 0}
};

void R_init_scInt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
