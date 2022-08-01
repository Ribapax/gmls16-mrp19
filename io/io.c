/*
 * Authors:
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *  Mateus Ribamar - GRR
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <dirent.h>
#include "../common/common.h"
#include "../linear_system/linear_system.h"

int ShowHelp(char *name) {
    fprintf(stderr, "\
        [usage] %s <options>\n\
        -e STRING   OPTIONAL  input file path. Default: stdin. \n\
        -s STRING   OPTIONAL  output file path. Default: stdout. \n\
        -r INTEGER  OPTIONAL  random matrix dimension. If present, ignores any input. \n\
        -i INTEGER  REQUIRED  refinement iterations limit (i >= 0).\n\
        ",
        name
    );
    return -1;
}

int GetOptions(
    int argc,
    char **argv,
    int *n,
    int *iterationsLimit,
    char **outputFilePath,
    char **inputFilePath
) {
    const struct option options[] = {
            {"e", optional_argument,  0, 'e'},
            {"s", optional_argument,  0, 's'},
            {"i", optional_argument,   0, 'i'},
            {"r", optional_argument,   0, 'r'},
            {0, 0, 0, 0},
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "he:s:i:r:", options, NULL)) > 0) {
        switch (opt) {
            case 'h':
                ShowHelp(argv[0]);
                break;
            case 'e':
                *inputFilePath = optarg;
                break;
            case 's':
                *outputFilePath = optarg;
                break;
            case 'i':
                *iterationsLimit = atoi(optarg);
                if (*iterationsLimit < 0) {
                    fprintf(stderr, "iterations limit should be greater or equal zero\n");
                    return -1;
                }
                break;
            case 'r':
                *n = atoi(optarg);
                break;
            default:
                return -1;
        }
    }
    return 0;
}

void PrintMatrix(FILE *outputFile, RealNumber **A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fprintf(outputFile, "%.15g ", A[i][j]);
        }
        fprintf(outputFile, "\n");
    }
}

RealNumber **ReadMatrix(char *fileName, int *size) {
    FILE *file = stdin;
    if (fileName) {
        file = fopen(fileName, "r");
        if (!file) {
            fprintf(stderr, "could not open file %s\n", fileName);
            return NULL;
        }
    }
    int n;
    fscanf(file, "%d", &n);
    *size = n;
    RealNumber **A = AllocateLinearSystem(n, PointerToPointer)->A;
    if (A == NULL) {
        fprintf(stderr, "could not allocate matrix\n");
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    fclose(file);
    return A;
}