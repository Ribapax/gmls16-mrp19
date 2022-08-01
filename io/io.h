/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR
 */

#ifndef IO_HEADER
#define IO_HEADER

#include "../common/common.h"
#include <stdio.h>

int ShowHelp(char *name);
int GetOptions(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath);
void PrintMatrix(FILE *outputFile, RealNumber **A, int n);
RealNumber **ReadMatrix(char *fileName, int *size);

#endif