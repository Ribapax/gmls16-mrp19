#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void showHelp(char *name) {
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

int getOptions(int argc, char **argv, int *n, int *iterationsLimit, char **outputFilePath, char **inputFilePath) {
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
                showHelp(argv[0]);
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