#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <dirent.h>
#include "../common/common.h"
#include "../linear_system/linear_system.h"

void ShowHelp(char *name) {
    fprintf(stderr, "\
        [usage] %s <options>\n\
        -e STRING   OPTIONAL  input file path.\n\
        -s STRING   OPTIONAL  output file.\n\
        -r INTEGER  OPTIONAL  matrix dimension.\n\
        -i INTEGER  REQUIRED  iterations limit (i > 0).\n\
        ",
        name
    );
    exit(-1);
}

int GetOptions(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath) {
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
                break;
            case 'r':
                *n = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Invalid option or missing argument: `%c'\n", optopt) ;
                return -1;
        }
    }
    return 0;
}

void printMatrix(RealNumber **A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%.15g ", A[i][j]);
        }
        printf("\n");
    }
}

RealNumber **ReadMatrix(char *fileName, int *size) {
    FILE *file = stdin;
    if (strlen(fileName) != 0) {
        file = fopen(fileName, "r");
    }
    if (!file) {
        perror("could not open file");
        exit(1);
    }
    int n;
    fscanf(file, "%d", &n);
    *size = n;
    RealNumber **A = AllocateLinearSystem(n, PointerToPointer)->A;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    fclose(file);
    return A;
}