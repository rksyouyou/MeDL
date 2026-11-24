#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP reml_full(SEXP, SEXP, SEXP);
extern SEXP lr_profiles(SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"reml_full",   (DL_FUNC) &reml_full,   3},
    {"lr_profiles", (DL_FUNC) &lr_profiles, 3},
    {NULL, NULL, 0}
};

void attribute_visible R_init_MeDL(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
