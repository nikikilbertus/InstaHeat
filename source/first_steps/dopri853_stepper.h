#ifndef __DOPRI853_STEPPER__
#define __DOPRI853_STEPPER__

#include "main.h"

typedef struct {
double c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c14,c15,c16,
b1,b6,b7,b8,b9,b10,b11,b12,
bhh1,bhh2,bhh3,
er1,er6,er7,er8,er9,er10,er11,er12,
a21,
a31,a32,
a41,a43,
a51,a53,a54,
a61,a64,a65,
a71,a74,a75,a76,
a81,a84,a85,a86,a87,
a91,a94,a95,a96,a97,a98,
a101,a104,a105,a106,a107,a108,a109,
a111,a114,a115,a116,a117,a118,a119,a1110,
a121,a124,a125,a126,a127,a128,a129,a1210,a1211,
a141,a147,a148,a149,a1410,a1411,a1412,a1413,
a151,a156,a157,a158,a1511,a1512,a1513,a1514,
a161,a166,a167,a168,a169,a1613,a1614,a1615,
d41,d46,d47,d48,d49,d410,d411,d412,d413,d414,d415,d416,
d51,d56,d57,d58,d59,d510,d511,d512,d513,d514,d515,d516,
d61,d66,d67,d68,d69,d610,d611,d612,d613,d614,d615,d616,
d71,d76,d77,d78,d79,d710,d711,d712,d713,d714,d715,d716;
}dopri853_constants_t;

typedef struct {
	double t;
    double t_old;
    double ti;
    double tf;
    double dt;
    double dt_did;
    double dt_next;
    double dt_min;
    size_t Ntot2;
    size_t max_steps;
    int n_stp;
    int n_ok;
    int n_bad;
    double beta, alpha;
    double safe;
    double minscale, maxscale;
	double a_tol, r_tol;
	double err_old;
	int reject;
	double eps;
    int dense;
}dopri853_control_t;

typedef struct {
	double *k2, *k3, *k4, *k5, *k6, *k7, *k8, *k9, *k10, *k_tmp;
	double *yerr, *yerr2;
	double *rcont1, *rcont2, *rcont3, *rcont4,
		   *rcont5, *rcont6, *rcont7, *rcont8;
	double a2, a3, a4, a5, a6, a7, a8, a9, a10, a_tmp;
}dopri853_values_t;

extern dopri853_constants_t dpc;
extern dopri853_control_t dp;
extern dopri853_values_t dpv;

void initialize_dopri853(parameters_t *pars);
void run_dopri853(parameters_t *pars);
void perform_step(const double dt_try, parameters_t *pars);
void try_step(const double dt, parameters_t *pars);
double error(const double dt);
int success(const double err, double *dt);
void prepare_dense_output(const double dt, parameters_t *pars);
double dense_output(const size_t i, const double t, const double dt);
void allocate_dopri853_values();
void destroy_dopri853_values();

#endif