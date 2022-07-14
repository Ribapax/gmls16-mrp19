#include <stdio.h>
#include <stdlib.h>
#include "io/io.h"

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

    printf("Hello, World!\n");
    return 0;
}