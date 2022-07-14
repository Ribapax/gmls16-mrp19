#include <stdio.h>
#include <string.h>

#include "common.h"

/* Returns the elapsed time in milliseconds
 * Usage:
 * double time;
 * time = timestamp();
 * <code snippet>
 * time = timestamp() - time;
*/
double timestamp(void) {
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    return((double)(tp.tv_sec*1.0e3 + tp.tv_nsec*1.0e-6));
}

/* Generates the string '<baseName>_n'
 * i.e. if baseName = "ABC" and n = 10,
 * it returns "ABC_10"
 * Useful for generating makers for LIKWID.
*/
String markerName(String baseName, int n) {
    String mark = (String) malloc((strlen(baseName) + 1) + numDigits(n) + 1 );
    sprintf(mark, "%s_%u", baseName,n);
    // printf("*** %s\n", mark);
    return mark;
}
