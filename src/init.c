#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _RHD_f_nd(SEXP);
extern SEXP _RHD_f_nd_each(SEXP);
extern SEXP _RHD_Fadj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_Fadj_inner(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_Fadj_inner_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_Fadj_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_FPCA(SEXP, SEXP, SEXP);
extern SEXP _RHD_IQRout(SEXP, SEXP);
extern SEXP _RHD_IQRout_vec(SEXP, SEXP);
extern SEXP _RHD_MC1nd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_MC1nd_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_MC2one(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_MC2one_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_rank_rcpp(SEXP, SEXP);
extern SEXP _RHD_RHD_inner(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_RHD_inner_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_RHDonly(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_RHDonly_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_RKHS(SEXP, SEXP, SEXP);
extern SEXP _RHD_sort_cpp(SEXP);
extern SEXP _RHD_test(SEXP);
extern SEXP _RHD_Xout(SEXP, SEXP, SEXP, SEXP);
extern SEXP _RHD_Xout_all(SEXP, SEXP, SEXP);
extern SEXP _RHD_Xsim(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_RHD_f_nd",            (DL_FUNC) &_RHD_f_nd,             1},
  {"_RHD_f_nd_each",       (DL_FUNC) &_RHD_f_nd_each,        1},
  {"_RHD_Fadj",            (DL_FUNC) &_RHD_Fadj,            11},
  {"_RHD_Fadj_inner",      (DL_FUNC) &_RHD_Fadj_inner,       8},
  {"_RHD_Fadj_inner_prob", (DL_FUNC) &_RHD_Fadj_inner_prob,  7},
  {"_RHD_Fadj_prob",       (DL_FUNC) &_RHD_Fadj_prob,        9},
  {"_RHD_FPCA",            (DL_FUNC) &_RHD_FPCA,             3},
  {"_RHD_IQRout",          (DL_FUNC) &_RHD_IQRout,           2},
  {"_RHD_IQRout_vec",      (DL_FUNC) &_RHD_IQRout_vec,       2},
  {"_RHD_MC1nd",           (DL_FUNC) &_RHD_MC1nd,           10},
  {"_RHD_MC1nd_prob",      (DL_FUNC) &_RHD_MC1nd_prob,      10},
  {"_RHD_MC2one",          (DL_FUNC) &_RHD_MC2one,          10},
  {"_RHD_MC2one_prob",     (DL_FUNC) &_RHD_MC2one_prob,     10},
  {"_RHD_rank_rcpp",       (DL_FUNC) &_RHD_rank_rcpp,        2},
  {"_RHD_RHD_inner",       (DL_FUNC) &_RHD_RHD_inner,       14},
  {"_RHD_RHD_inner_prob",  (DL_FUNC) &_RHD_RHD_inner_prob,  13},
  {"_RHD_RHDonly",         (DL_FUNC) &_RHD_RHDonly,         12},
  {"_RHD_RHDonly_prob",    (DL_FUNC) &_RHD_RHDonly_prob,    11},
  {"_RHD_RKHS",            (DL_FUNC) &_RHD_RKHS,             3},
  {"_RHD_sort_cpp",        (DL_FUNC) &_RHD_sort_cpp,         1},
  {"_RHD_test",            (DL_FUNC) &_RHD_test,             1},
  {"_RHD_Xout",            (DL_FUNC) &_RHD_Xout,             4},
  {"_RHD_Xout_all",        (DL_FUNC) &_RHD_Xout_all,         3},
  {"_RHD_Xsim",            (DL_FUNC) &_RHD_Xsim,             5},
  {NULL, NULL, 0}
};

void R_init_RHD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
