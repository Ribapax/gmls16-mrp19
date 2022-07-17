#ifndef IO_HEADER
#define IO_HEADER

#include "../common/common.h"

void ShowHelp(char *name);
int GetOptions(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath);
void printMatrix(RealNumber **A, int n);
RealNumber **ReadMatrix(char *fileName, int *size);

#endif