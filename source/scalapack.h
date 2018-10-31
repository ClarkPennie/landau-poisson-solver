#ifndef LAPACK_H
#define LAPACK_H

extern "C" void dgetrf(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
extern "C" void dgetri(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);

#endif
