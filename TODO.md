# TODO

- [x] Develop and test `InvertMatrix` function (remember to remove the stub matrix from there);
- [x] Return errors as integer codes, and fail gracefully logging the error in the `stderr`;
- [x] Return error when the matrix cannot be inverted (determinant = 0). After calculating `U`, we need
    to test if the determinant is zero; just multiply all elements of the main diagonal of the matrix; if its zero, return an error.
- [x] Refinement.
- [ ] Write `LEIAME` file;
- [ ] Fix stopping criteria;
- [ ] Test the code against giant matrices;
- [ ] Test the code against non-reversible matrices;
- [ ] Nice to have: we can create enums for all the errors and write a function
called `handleError(int err)` that reads the enum, print the correct message to `stderr`
and `exit`;

## IMPORTANT
- [ ] The solution is changing when the turn on the Partial Pivoting. It should not change;
- [ ] Residue is not getting smaller with the iterations;
- [ ] Program does not return the correct result when doing more than 0 iterations.