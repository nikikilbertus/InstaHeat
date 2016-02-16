#ifndef __FULL_MULTIGRID__
#define __FULL_MULTIGRID__

#include <stddef.h>

#define NPRE            (1)
#define NPOST           (1)
#define NGMAX           (15)

void mglin(double **u, size_t n, size_t ncycle);
void addint(double **uf, double **uc, double **res, size_t nf);
void copy(double **aout, double **ain, size_t n);
void fill0(double **u, size_t n);
void interp(double **uf, double **uc, size_t nf);
void relax(double **u, double **rhs, size_t n);
void resid(double **res, double **u, double **rhs, size_t n);
void rstrct(double **uc, double **uf, size_t nc);
void slvsml(double **u, double **rhs);
void malloc_dmat(double*** mat, size_t n, size_t m);
void free_dmat(double*** mat, size_t n);
#endif
