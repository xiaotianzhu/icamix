#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP icamix_EMInterwovenFastICA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP icamix_WtsFastICA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP icamix_WtsKde(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"icamix_EMInterwovenFastICA", (DL_FUNC) &icamix_EMInterwovenFastICA, 8},
    {"icamix_WtsFastICA",          (DL_FUNC) &icamix_WtsFastICA,          7},
    {"icamix_WtsKde",              (DL_FUNC) &icamix_WtsKde,              4},
    {NULL, NULL, 0}
};

void R_init_icamix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
