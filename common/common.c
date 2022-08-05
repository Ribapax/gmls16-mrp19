/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <string.h>

#include "common.h"

double GetTimestamp(void) {
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    return((double)(tp.tv_sec*1.0e3 + tp.tv_nsec*1.0e-6));
}

String markerName(String baseName, int n) {
    String mark = (String) malloc((strlen(baseName) + 1) + numDigits(n) + 1 );
    sprintf(mark, "%s_%u", baseName,n);
    // printf("*** %s\n", mark);
    return mark;
}
