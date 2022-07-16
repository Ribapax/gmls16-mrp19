#include <stdio.h>
#include <stdlib.h>
#include "io/io.h"
#include "common/common.h"
#include "linear_system/linear_system.h"
#include "linear_system/lu_factorization.h"

int main(int argc, char *argv[]) {
    srand(20221);
    if (argc < 2) {
        showHelp(argv[0]);
    }

    int n;
    int iterationsLimit;
    char *outputFilePath, *inputFilePath;
    getOptions(argc, argv, &n, &iterationsLimit, &outputFilePath, &inputFilePath);

    printf("Params: \n");
    printf("n: %d\n", n);
    printf("iterationsLimit: %d\n", iterationsLimit);
    printf("outputFilePath: %s\n", outputFilePath);
    printf("inputFilePath: %s\n", inputFilePath);

    // #### Generate L matrix test ####

    int size = 3;

    // Matrix A to be inverted
    RealNumber **A = allocateLinearSystem(size, PointerToPointer)->A;
    A[0][0] = 25;
    A[0][1] = 5;
    A[0][2] = 1;
    A[1][0] = 64;
    A[1][1] = 8;
    A[1][2] = 1;
    A[2][0] = 144;
    A[2][1] = 12;
    A[2][2] = 1;

    RealNumber **U = allocateLinearSystem(size, PointerToPointer)->A;
    RealNumber **B = getIdentityMatrix(size);
    RealNumber **L = generateMatrixL(A, B, U, size);

    // Printing the multiplier's matrix L
    printf("\n L matrix: \n");
    printMatrix(L, size);
    printf("\n");

    // #### Final Code Algorithm ####

    // Read input

    // For each matrix
    // 1) invert the matrix
    // 2) calculate and print Residue L2 Norm
    // 3) print times
    // 4) print size and inverted matrix

    return 0;
}