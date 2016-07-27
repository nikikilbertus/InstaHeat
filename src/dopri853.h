#ifndef __DOPRI853_STEPPER__
#define __DOPRI853_STEPPER__

#include "main.h"
#if INTEGRATION_METHOD == DOPRI853

/**
 * @file dopri853.h
 *
 * @brief Struct definitions and function declarations for `dopri853.c` and for
 * `dopri853_constants.c`.
 *
 * @see <a href="http://numerical.recipes">Numerical Recipes</a>
 */

/**
 * @brief The Butcher tableaux for the Dormand Prince 8(5,3) (DOPRI853)
 * integration routine.
 */
struct dopri853_constants
{
    const double c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c14,c15,c16,
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
};

/**
 * @brief The only instance of `dopri853_constants`.
 */
extern struct dopri853_constants dpc;

void run_dopri853();

#endif
#endif
