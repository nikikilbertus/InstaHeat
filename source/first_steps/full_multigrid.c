#include <stdlib.h>
#include "full_multigrid.h"
#include "main.h"

void mglin(double **u, size_t n, size_t ncycle) {
    size_t j, jcycle, jj, jpost, jpre, nf, ng = 0, ngrid, nn;
    double **ires[NGMAX + 1], **irho[NGMAX + 1], **irhs[NGMAX + 1], **iu[NGMAX + 1];

    nn = n;
    while (nn >>= 1)
    {
        ng += 1;
    }

    if (n != 1 + (1L << ng))
    {
        fputs("n-1 must be a power of 2 in full multigrid", stderr);
        exit(EXIT_FAILURE);
    }

    if (ng > NGMAX)
    {
        fputs("increase NGMAX for full multigrid", stderr);
        exit(EXIT_FAILURE);
    }

    nn = n / 2 + 1;
    ngrid = ng - 1;
    malloc_dmat(&(irho[ngrid]), nn, nn);
    rstrct(irho[ngrid], u, nn);
    while (nn > 3)
    {
        nn = nn / 2 + 1;
        malloc_dmat(&(irho[--ngrid]), nn, nn);
        rstrct(irho[ngrid], irho[ngrid + 1], nn);
    }

    nn = 3;
    malloc_dmat(&(iu[1]), nn, nn);
    malloc_dmat(&(irhs[1]), nn, nn);
    slvsml(iu[1], irho[1]);
    free_dmat(&(irho[1]), nn);
    ngrid = ng;
    for (j = 2; j <= ngrid; ++j)
    {
        nn = 2 * nn - 1;
        malloc_dmat(&(iu[j]), nn, nn);
        malloc_dmat(&(irhs[j]), nn, nn);
        malloc_dmat(&(ires[j]), nn, nn);
        interp(iu[j], iu[j - 1], nn);
        copy(irhs[j], (j != ngrid ? irho[j] : u), nn);

        for (jcycle = 1; jcycle <= ncycle; ++jcycle)
        {
            nf = nn;
            for (jj = j; jj >= 2; --jj)
            {
                for (jpre = 1; jpre <= NPRE; ++jpre)
                {
                    relax(iu[jj], irhs[jj], nf);
                }
                resid(ires[jj], iu[jj], irhs[jj], nf);
                nf = nf / 2 + 1;
                rstrct(irhs[jj - 1], ires[jj], nf);
                fill0(iu[jj - 1], nf);
            }
            slvsml(iu[1], irhs[1]);
            nf = 3;
            for (jj = 2; jj <= j; ++jj)
            {
                nf = 2 * nf - 1;
                addint(iu[jj], iu[jj - 1], ires[jj], nf);
                for (jpost = 1; jpost <= NPOST; ++jpost)
                {
                    relax(iu[jj], irhs[jj], nf);
                }
            }
        }
    }
    copy(u, iu[ngrid], n);
    for (nn = n, j = ng; j >= 2; --j, nn = nn / 2 + 1)
    {
        free_dmat(&(ires[j]), nn);
        free_dmat(&(irhs[j]), nn);
        free_dmat(&(iu[j]), nn);
        if (j != ng)
        {
            free_dmat(&(irho[j]), nn);
        }
    }
    free_dmat(&(irhs[1]), 3);
    free_dmat(&(iu[1]), 3);
}

void rstrct(double **uc, double **uf, size_t nc) {
    size_t ic, iif, jc, jf, ncc = 2 * nc - 1;

    for (jf = 3, jc = 2; jc < nc; ++jc, jf += 2)
    {
        for (iif = 3, ic = 2; ic < nc; ++ic, iif += 2)
        {
            uc[ic][jc] = 0.5 * uf[iif][jf] + 0.125 * (uf[iif + 1][jf] +
                    uf[iif - 1][jf] + uf[iif][jf + 1] + uf[iif][jf - 1]);
        }
    }

    for (jc = 1, ic = 1; ic <= nc; ++ic, jc += 2)
    {
        uc[ic][1]  = uf[jc][1];
        uc[ic][nc] = uf[jc][ncc];
    }

    for (jc = 1, ic = 1; ic <= nc; ++ic, jc += 2)
    {
        uc[1][ic]  = uf[1][jc];
        uc[nc][ic] = uf[ncc][ic];
    }
}

void interp(double **uf, double **uc, size_t nf) {
    size_t ic, iif, jc, jf, nc;
    nc = nf / 2 + 1;
    for (jc = 1, jf = 1; jc <= nc; ++jc, jf += 2)
    {
        for (ic = 1; ic <= nc; ++ic)
        {
            uf[2 * ic - 1][jf] = uc[ic][jc];
        }
    }

    for (jf = 1; jf <= nf; jf += 2)
    {
        for (iif = 2; iif < nf; iif += 2)
        {
            uf[iif][jf] = 0.5 * (uf[iif][jf+1] + uf[iif][jf-1]);
        }
    }

    for (jf = 2; jf < nf; jf += 2)
    {
        for (iif = 1; iif <= nf; ++iif)
        {
            uf[iif][jf] = 0.5 * (uf[iif][jf+1] + uf[iif][jf-1]);
        }
    }
}

void addint(double **uf, double **uc, double **res, size_t nf) {
    size_t i, j;
    interp(res, uc, nf);
    for (j = 1; j <= nf; ++j)
    {
        for (i = 1; i <= nf; ++i)
        {
            uf[i][j] += res[i][j];
        }
    }
}

//TODO: modify for our problem after testing
void slvsml(double **u, double **rhs) {
    //TODO: change!
    double h = 0.5;

    fill0(u,3);
    u[2][2] = - h * h * rhs[2][2] / 4.0;
}

void relax(double **u, double **rhs, size_t n) {
    size_t i, ipass, isw, j, jsw = 1;
    double h, h2;

    //TODO change
    h = 1.0 / (n-1);
    h2 = h * h;

    for (ipass = 1; ipass <= 2; ++ipass, jsw = 3 - jsw)
    {
        isw = jsw;
        for (j = 2; j < n; ++j, isw = 3 - isw)
        {
            for (i = isw + 1; i < n; i += 2)
            {
                u[i][j] = 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] +
                        u[i][j - 1] - h2 * rhs[i][j]);
            }
        }
    }
}

void resid(double **res, double **u, double **rhs, size_t n) {
    size_t i, j;
    double h, h2i;

    //TODO: change
    h = 1.0 / (n - 1);
    h2i = 1.0 / (h * h);

    for (j = 2; j < n; ++j)
    {
        for (i = 2; i < n; ++i)
        {
            res[i][j] = -h2i * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] +
                    u[i][j - 1] - 4.0 * u[i][j]) + rhs[i][j];
        }
    }

    for (i = 1; i <= n; ++i)
    {
        res[i][1] = res[i][n] = res[1][i] = res[n][i] = 0.0;
    }
}

void copy(double **aout, double **ain, size_t n) {
    size_t i, j;
    for (i = 1; i <= n; ++i)
    {
        for (j = 1; j <= n; ++j)
        {
            aout[j][i] = ain[j][i];
        }
    }
}

void fill0(double **u, size_t n) {
    size_t i, j;
    for (j = 1; j <= n; ++j)
    {
        for (i = 1; i <= n; ++i)
        {
            u[i][j] = 0.0;
        }
    }
}

void malloc_dmat(double*** mat, size_t n, size_t m) {
    *mat = (double**) malloc(n * sizeof(double*));
    for (size_t i = 0; i < n; ++i)
    {
        (*mat)[i] = (double*) malloc(m * sizeof(double));
    }
}

void free_dmat(double*** mat, size_t n) {
    for (size_t i = 0; i < n; ++i)
    {
        free((*mat)[i]);
    }
    free(*mat);
}
