/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR20190154
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/io.h"
#include "common/common.h"
#include "linear_system/linear_system.h"
#include "linear_system/lu_factorization.h"


int main(int argc, char *argv[]) {
    srand(S_RAND_CONST);
    int err;
    if (argc < 2) {
        err = ShowHelp(argv[0]);
        if (err != 0) {
            exit(-1);
        }
    }
    int randomMatrixSize = 0;
    int iterationsLimit;
    char *outputFilePath = NULL, *inputFilePath = NULL;
    err = GetOptions(
        argc,
        argv,
        &randomMatrixSize,
        &iterationsLimit,
        &outputFilePath,
        &inputFilePath
    );
    if (err != 0) {
        fprintf(stderr, "invalid option or missing argument\n");
        exit(-1);
    }

    FILE *outputFile = stdout;
    if (outputFilePath) {
        outputFile = fopen(outputFilePath, "w+");
        if (!outputFile) {
            fprintf(stderr, "could not open file %s", outputFilePath);
            exit(-1);
        }
        printf("\nResults stored in the file \"%s\".\n", outputFilePath);
    }

    // Read input
    int size;
    RealNumber **A = NULL;
    if (randomMatrixSize) {
        LinearSystem *LS = AllocateLinearSystem(randomMatrixSize, PointerToPointer);
        if (LS == NULL) {
            fprintf(stderr, "could not allocate matrix\n");
            exit(-1);
        }
        err = FillLinearSystem(LS, GenericMatrix, COEFFICIENT_LIMIT);
        if (err != 0) {
            fprintf(stderr, "could not fill matrix with random values\n");
            exit(-1);
        }
        A = LS->A;
        size = randomMatrixSize;
    } else {
        A = ReadMatrix(inputFilePath, &size);
        if (A == NULL) {
            fprintf(stderr, "could not read matrix\n");
            exit(-1);
        }
    }

    if (!MatrixIsInvertible(A, size)) {
        fprintf(stderr, "matrix is not invertible");
        exit(-1);
    }

    fprintf(outputFile, "#\n");

    RealNumber **B = GetIdentityMatrix(size);
    if (B == NULL) {
        fprintf(stderr, "could not allocate identity matrix B\n");
        exit(-1);
    }

    RealNumber **U = AllocateLinearSystem(size, PointerToPointer)->A;
    if (U == NULL) {
        fprintf(stderr, "could not allocate \"U\" matrix\n");
        exit(-1);
    }

    RealNumber **L = AllocateLinearSystem(size, PointerToPointer)->A;
    if (L == NULL) {
        fprintf(stderr, "could not allocate \"L\" matrix\n");
        exit(-1);
    }

    Time avgLSTime = 0, LUTime = 0, residueTime = 0;
    RealNumber **invertedA = InvertMatrix(A, B, size, &avgLSTime, &LUTime, L, U);
    if (invertedA == NULL) {
        fprintf(stderr, "could not invert matrix\n");
        exit(-1);
    }

    Time residueTimeTemp = timestamp();
    RealNumber currentResidueL2Norm = CalculateResidue(A, B, invertedA, size);
    residueTime += timestamp() - residueTimeTemp;
    fprintf(outputFile, "# iter 1: %.15g\n", currentResidueL2Norm);

    RealNumber lastResidueL2Norm;
    int iteration = 2;
    do {
        lastResidueL2Norm = currentResidueL2Norm;
        invertedA = refineSolution(A, B, invertedA, L, U, size);
        residueTimeTemp = timestamp();
        currentResidueL2Norm = CalculateResidue(A, B, invertedA, size);
        residueTime += timestamp() - residueTimeTemp;
        if (ResidueIsIncreasing(currentResidueL2Norm, lastResidueL2Norm)) {
            continue;
        }
        fprintf(outputFile, "# iter %d: %.15g\n", iteration, currentResidueL2Norm);
        iteration++;
    } while (
        hasNotReachedStoppingCriteria(iteration, iterationsLimit, currentResidueL2Norm, lastResidueL2Norm)
    );

    fprintf(outputFile, "# Tempo LU: %.15g\n", LUTime);
    fprintf(outputFile, "# Tempo iter: %.15g\n", avgLSTime);
    fprintf(outputFile, "# Tempo residuo: %.15g\n", residueTime);
    fprintf(outputFile, "#\n");
    fprintf(outputFile, "%d\n", size);
    PrintMatrix(outputFile, invertedA, size);

    return 0;
}