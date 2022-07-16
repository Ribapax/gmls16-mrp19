#ifndef IO_HEADER
#define IO_HEADER

#include "../common/common.h"

void showHelp(char *name);
int getOptions(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath);
void printMatrix(RealNumber **A, int n);

#endif