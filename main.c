#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/io.h"
#include "common/common.h"
#include "linear_system/linear_system.h"
#include "linear_system/lu_factorization.h"

RealNumber **refineSolution(RealNumber **A, RealNumber **invertedMatrix, int n) {
    RealNumber *Ax;
    RealNumber *B;
    RealNumber **refinedSolution = AllocateLinearSystem(n, PointerToPointer)->A;
    RealNumber **identityMatrix = GetIdentityMatrix(n);
    RealNumber *currentColumn = malloc(sizeof(RealNumber) * n);
    RealNumber *currentIdentityColumn = malloc(sizeof(RealNumber) * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; ++j) {
            currentColumn[j] = invertedMatrix[j][i];
            currentIdentityColumn[j] = identityMatrix[j][i];
        }
        Ax = multiplyMatrixWithArray(A, currentColumn, n);
        B = subtractArrays(currentIdentityColumn, Ax, n);
        refinedSolution[i] = GaussElimination(A, B, n);
    }
    return refinedSolution;
}

int main(int argc, char *argv[]) {
    srand(S_RAND_CONST);
    int err = 0;
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
    RealNumber **invertedA;
    Time residueTime = 0.;
    Time residueTimeTemp;
    Time avgLSTime = 0, LUTime = 0;
    RealNumber residueL2Norm = 1 + RESIDUE_THRESHOLD;

    invertedA = InvertMatrix(A, size, &avgLSTime, &LUTime);
    if (invertedA == NULL) {
        fprintf(stderr, "could not invert matrix\n");
        exit(-1);
    }
    residueTimeTemp = timestamp();
    residueL2Norm = CalculateResidue(A, invertedA, size);
    residueTime += timestamp() - residueTimeTemp;
    int iteration = 1;
    while (hasNotReachedStoppingCriteria(iteration, iterationsLimit, residueL2Norm)) {
    // 1) Refine Solution
        invertedA = refineSolution(A, invertedA, size);
    // 2) Calculate Residue L2 Norm
        residueTimeTemp = timestamp();
        residueL2Norm = CalculateResidue(A, invertedA, size);
        residueTime += timestamp() - residueTimeTemp;
        fprintf(outputFile, "# iter %d: %.15g\n", iteration, residueL2Norm);
        iteration++;
    }

    // 3) print times
    fprintf(outputFile, "# Tempo LU: %.15g\n", LUTime);
    fprintf(outputFile, "# Tempo iter: %.15g\n", avgLSTime);
    fprintf(outputFile, "# Tempo residuo: %.15g\n", residueTime);
    fprintf(outputFile, "#\n");
    // 4) print size and inverted matrix
    fprintf(outputFile, "%d\n", size);
    PrintMatrix(outputFile, invertedA, size);

    return 0;
}