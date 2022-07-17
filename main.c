#include <stdio.h>
#include <stdlib.h>
#include "io/io.h"
#include "common/common.h"
#include "linear_system/linear_system.h"
#include "linear_system/lu_factorization.h"

int main(int argc, char *argv[]) {
    srand(S_RAND_CONST);
    if (argc < 2) {
        ShowHelp(argv[0]);
    }
    int randomMatrixSize = 0;
    int iterationsLimit;
    char *outputFilePath, *inputFilePath;
    GetOptions(argc, argv, &randomMatrixSize, &iterationsLimit, &outputFilePath, &inputFilePath);

    // Read input
    int size;
    RealNumber **A;
    if (randomMatrixSize) {
        LinearSystem *LS = AllocateLinearSystem(randomMatrixSize, PointerToPointer);
        FillLinearSystem(LS, GenericMatrix, COEFFICIENT_LIMIT);
        A = LS->A;
        size = randomMatrixSize;
    } else {
        A = ReadMatrix(inputFilePath, &size);
    }

    printf("#\n");
    RealNumber **invertedA;
    Time residueTime;
    Time avgLSTime = 0, LUTime = 0;
    RealNumber residueL2Norm = 9999;
    int iteration = 1;
    while (hasNotReachedStoppingCriteria(iteration, iterationsLimit, residueL2Norm)) {
        // 1) invert the matrix
        invertedA = InvertMatrix(A, size, &avgLSTime, &LUTime);
        // 2) calculate and print Residue L2 Norm
        residueTime = timestamp();
        residueL2Norm = CalculateResidue(A, invertedA, size);
        residueTime = timestamp() - residueTime;
        printf("# iter %d: %.15g\n", iteration, residueL2Norm);
        iteration++;
    }

    // 3) print times
    printf("# Tempo LU: %.15g\n", LUTime);
    printf("# Tempo iter: %.15g\n", avgLSTime);
    printf("# Tempo residuo: %.15g\n", residueTime);
    printf("#\n");
    // 4) print size and inverted matrix
    printf("%d\n", size);
    printMatrix(invertedA, size);

    return 0;
}