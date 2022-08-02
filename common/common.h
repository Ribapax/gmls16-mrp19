/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef COMMON_HEADER
#define COMMON_HEADER

#include <stdlib.h>
#include <time.h>

#define S_RAND_CONST 20221

// Absolute value of a number, alternative solution to fabs()
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

typedef double RealNumber;
typedef char * String;
typedef double Time;

// SIMD alignment macros
#define ALIGN_64 __attribute__((aligned(64)))
#define ALIGN_32 __attribute__((aligned(32)))
#define ALIGN_16 __attribute__((aligned(16)))

// Limit of digits in a number
#define numDigits(n)  6  // ( (int) log10(n) + 1 )

// Verify if n is power of 2
#define isPot2(n) ((n) && !((n) & ((n) - 1)))

/* Returns the elapsed time in milliseconds
 * Usage:
 * double time;
 * time = timestamp();
 * <code snippet>
 * time = timestamp() - time;
*/
double timestamp(void);

/* Generates the string '<baseName>_n'
 * i.e. if baseName = "ABC" and n = 10,
 * it returns "ABC_10"
 * Useful for generating makers for LIKWID.
*/
String markerName(String baseName, int n);

#endif // COMMON_HEADER

