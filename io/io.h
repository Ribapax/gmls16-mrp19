/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#ifndef IO_HEADER
#define IO_HEADER

#include "../common/common.h"
#include <stdio.h>

/*
 * Prints a usage guide to stdout and returns an error value.
 * Params:
 *  - name: the name of the program;
 * Always Returns -1 (error).
 */
int ShowHelp(char *name);

/*
 * Get the options passed as argument to the program.
 * Params:
 *  - n: random matrix dimension;
 *  - iterationsLimit: maximum number of iterations of the refinement method;
 *  - outputFilePath: path of the file in which the results will be stored;
 *  - inputFilePath: path to the file containing the input matrix.
 * Returns:
 *  - 0 (success);
 *  - -1 (error).
 */
int ReadParameters(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath);

/*
 * Prints a matrix `A` of size `n` to `outputFile`.
 * Params:
 *  - outputFile: file path or "stdout";
 *  - A: matrix;
 *  - n: matrix dimension.
 * Returns nothing.
 */
void PrintMatrix(FILE *outputFile, RealNumber *A, int n);

/*
 * Reads a matrix of size `n` coming from `fileName`.
 * Params:
 *  - fileName: input file path or "stdin";
 *  - size: matrix dimension.
 * Returns:
 *  - NULL (some error occurred);
 *  - RealNumber** (a matrix correctly allocated).
 */
RealNumber *ReadMatrix(char *fileName, int *size);

#endif