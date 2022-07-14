#ifndef COMMON_HEADER
#define COMMON_HEADER

#include <stdlib.h>
#include <time.h>

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

double timestamp(void);
String markerName(String baseName, int n);

#endif // COMMON_HEADER

