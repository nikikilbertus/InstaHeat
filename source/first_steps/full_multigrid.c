#include "main.h"

void mglin(double **u, size_t n, size_t ncycle) {
    size_t j, jcycle, jj, jpost, jpre, nf, ng = 0, ngrid, nn;
    double **ired[NGMAX + 1], **irho[NGMAX + 1], **iu[NGMAX + 1];

    nn = n;

    while (nn >>= 1)
    {
        ng += 1;
    }

    if (n != 1 + (1L << ng))
    {
        //TODO: exit with error
    }

    if (ng > NGMAX)
    {
        //TODO: exit with error
    }

    nn = n/2 + 1;
    ngrid = ng - 1;
    irho[ngrid] = malloc(nn * nn * sizeof *irho[ngrid]);
    rstrct(irho[ngrid], u, nn);
    while (nn > 3)
    {
        nn = nn/2 + 1;
        irho[--ngrid] = malloc(nn * nn * sizeof *irho[ngrid]);
        rstrct(irho[ngrid], irho[ngrid + 1], nn);
    }

    nn = 3;
    iu[1]   = malloc(nn * nn * sizeof *iu[1]);
    irhs[1] = malloc(nn * nn * sizeof *irhs[1]);
    slvsml(iu[1], irho[1]);
    free(irho[1]);
    ngrid = ng;
    for (j = 2; j <= ngrid; ++j)
    {
        nn = 2 * nn - 1;
        iu[j]   = malloc(nn * nn * sizeof *iu[j]);
        irhs[j] = malloc(nn * nn * sizeof *irhs[j]);
        ires[j] = malloc(nn * nn * sizeof *ires[j]);
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
        free(ires[j]);
        free(irhs[j]);
        free(iu[j]);
        if (j != ng)
        {
            free(irho[j]);
        }
    }
    free(irhs[1]);
    free(iu[1]);
}
