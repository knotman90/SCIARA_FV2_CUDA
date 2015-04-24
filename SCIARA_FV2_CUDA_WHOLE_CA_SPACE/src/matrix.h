#ifndef matrix_h
#define matrix_h

#include <stdio.h>

//---------------------------------------------------------------------------
//MEMORY ALLOCATION
//---------------------------------------------------------------------------

template <class Type>
Type** allocateMatrix(Type** M, int m, int n);

template <class Type>
Type*** allocateMatrix(Type*** M, int m, int n, int k);

template <class Type>
void deAllocateMatrix(Type** M, int m, int n);

template <class Type>
void deAllocateMatrix(Type*** M, int m, int n, int k);

//---------------------------------------------------------------------------
//SET TO ZERO
//---------------------------------------------------------------------------

template <class Type>
Type** setToZeroMatrix(Type** M, int lx, int ly);

template <class Type>
Type*** setToZeroMatrix(Type*** M, int lx, int ly, int k);

template <class Type>
Type** setToZeroMatrix(Type** M, int lx_min, int lx_max, int ly_min, int ly_max);

template <class Type>
Type*** setToZeroMatrix(Type*** M, int lx_min, int lx_max, int ly_min, int ly_max, int k);

//---------------------------------------------------------------------------
//COPY
//---------------------------------------------------------------------------

template <class Type>
void copyMatrix(Type** M_Destinatario, Type** M_Sorgente, int lx, int ly);

template <class Type>
void copyMatrix(Type** M_Destinatario, Type** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max);

//---------------------------------------------------------------------------
//SUBTRACTION
//---------------------------------------------------------------------------

template <class Type>
void addMatrices(Type** Dest, Type** M1, Type** M2, int lx, int ly);

template <class Type>
void subtractMatrix(Type** Dest, Type** M1, Type** nenoM2, int lx, int ly);

//---------------------------------------------------------------------------
//INPUT - OUTPUT
//---------------------------------------------------------------------------

bool**		readMatrix(bool** M, int lx, int ly, FILE* f);
int**		readMatrix(int** M, int lx, int ly, FILE* f);
float**		readMatrix(float** M, int lx, int ly, FILE* f);
double**	readMatrix(double** M, int lx, int ly, FILE* f);

void saveMatrix(bool** M, int lx, int ly, FILE* f);
void saveMatrix(int** M, int lx, int ly, FILE* f);
void saveMatrix(float** M, int lx, int ly, FILE* f);
void saveMatrix(double** M, int lx, int ly, FILE* f);
void saveMatrix(int** M, int lx, int ly, int threshold, FILE* f);
void saveMatrix(double** M, int lx, int ly, double threshold, FILE* f);

template <class Type> Type** readMatrix(Type** M, int lx, int ly, char* path);
template <class Type> void saveMatrix(Type** M, int lx, int ly, char* path);
template <class Type> void saveMatrix(Type** M, int lx, int ly, Type threshold, char* path);

#endif

