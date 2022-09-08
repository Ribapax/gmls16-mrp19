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
#include "lu_factorization/lu_factorization.h"

// This block enables to compile the code with and without the LIKWID header in place
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

int main(int argc, char *argv[]) {

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;

    srand(S_RAND_CONST);
    int err;
    if (argc < 2) {
        err = ShowHelp(argv[0]);
        if (err != 0) {
            exit(-1);
        }
    }
    // Read parameters
    int randomMatrixSize = 0;
    int iterationsLimit;
    char *outputFilePath = NULL, *inputFilePath = NULL;
    err = ReadParameters(
        argc,
        argv,
        &randomMatrixSize,
        &iterationsLimit,
        &outputFilePath,
        &inputFilePath
    );
    if (err != 0) {
        fprintf(stderr, "invalid parameter or missing argument\n");
        exit(-1);
    }

    // Defining where to store or print the results
    FILE *outputFile = stdout;
    if (outputFilePath) {
        outputFile = fopen(outputFilePath, "w+");
        if (!outputFile) {
            fprintf(stderr, "could not open file %s", outputFilePath);
            exit(-1);
        }
        printf("\nResults stored in the file \"%s\".\n", outputFilePath);
    }

    // Reading the input Matrix
    int size;
    RealNumber **A = NULL;
    if (randomMatrixSize) {
        // Generate a Random Matrix
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
        // Read the Matrix from `stdin` or from the given file path
        A = ReadMatrix(inputFilePath, &size);
        if (A == NULL) {
            fprintf(stderr, "could not read matrix\n");
            exit(-1);
        }
    }

    

    fprintf(outputFile, "#\n");

    // Generate the Identity Matrix
    RealNumber **B = GenerateIdentityMatrix(size);
    if (B == NULL) {
        fprintf(stderr, "could not allocate identity matrix B\n");
        exit(-1);
    }

    // Allocate the Upper Matrix
    RealNumber **U = AllocateLinearSystem(size, PointerToPointer)->A;
    if (U == NULL) {
        fprintf(stderr, "could not allocate \"U\" matrix\n");
        exit(-1);
    }

    // Allocate the Lower Matrix
    RealNumber **L = AllocateLinearSystem(size, PointerToPointer)->A;
    if (L == NULL) {
        fprintf(stderr, "could not allocate \"L\" matrix\n");
        exit(-1);
    }

    // Invert the given Matrix
    Time avgLSTime = 0, LUTime = 0, residueTime = 0;
    LUTime = GetTimestamp();
    LUDecomposition(A, B, U, L, size);
    LUTime = GetTimestamp() - LUTime;

    // Test if the given Matrix is invertible
    if (!MatrixIsInvertible(U, size)) {
        fprintf(stderr, "matrix is not invertible");
        exit(-1);
    }

    LIKWID_MARKER_START("LINEAR_SYSTEM_CALCULATION");
    RealNumber **invertedA = SolveLinearSystems(B, size, &avgLSTime, L, U);
    LIKWID_MARKER_STOP("LINEAR_SYSTEM_CALCULATION");
    if (invertedA == NULL) {
        fprintf(stderr, "could not invert matrix\n");
        exit(-1);
    }

    // Calculate the first L2 Norm of the Residue
    Time residueTimeTemp = GetTimestamp();
    RealNumber currentResidueL2Norm = CalculateResidueL2Norm(A, B, invertedA, size);
    residueTime += GetTimestamp() - residueTimeTemp;
    fprintf(outputFile, "# iter 1: %.15g\n", currentResidueL2Norm);

    RealNumber lastResidueL2Norm;
    int iteration = 2;
    // Start refining the solution until it reaches the stopping criteria
    do {
        lastResidueL2Norm = currentResidueL2Norm;
        invertedA = RefineSolution(A, B, invertedA, L, U, &avgLSTime, size);
        residueTimeTemp = GetTimestamp();
        currentResidueL2Norm = CalculateResidueL2Norm(A, B, invertedA, size);
        residueTime += GetTimestamp() - residueTimeTemp;
        if (ResidueIsIncreasing(currentResidueL2Norm, lastResidueL2Norm)) {
            continue;
        }
        fprintf(outputFile, "# iter %d: %.15g\n", iteration, currentResidueL2Norm);
        iteration++;
    } while (
        HasNotReachedStoppingCriteria(iteration, iterationsLimit, currentResidueL2Norm, lastResidueL2Norm)
    );

    // Print the final timestamps
    fprintf(outputFile, "# Tempo LU: %.15g\n", LUTime);
    fprintf(outputFile, "# Tempo iter: %.15g\n", avgLSTime);
    fprintf(outputFile, "# Tempo residuo: %.15g\n", residueTime);
    fprintf(outputFile, "#\n");
    fprintf(outputFile, "%d\n", size);
    // Print the final Inverted Matrix
    PrintMatrix(outputFile, invertedA, size);

    LIKWID_MARKER_CLOSE;

    return 0;
}