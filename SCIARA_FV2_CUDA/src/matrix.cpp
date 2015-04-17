#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//---------------------------------------------------------------------------
//MEMORY ALLOCATION
//---------------------------------------------------------------------------

template <class Type>
Type** allocateMatrix(Type** M, int m, int n)
{
    M = new Type*[m];
    for (int i = 0; i < m; i++)
        M[i] = new Type[n];
    return M;
}

template <class Type>
Type*** allocateMatrix(Type*** M, int m, int n, int k)
{
    M = new Type**[m];
    for (int i = 0; i < m; i++)
    {
        M[i] = new Type*[n];
        for (int j = 0; j < n; j++)
            M[i][j] = new Type[k];
    }

    return M;
}

template <class Type>
void deAllocateMatrix(Type** M, int m, int n)
{
    for (int i = 0; i < m; i++)
        delete[] M[i];
    delete[] M;
    return;
}

template <class Type>
void deAllocateMatrix(Type*** M, int m, int n, int k)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            delete[] M[i][j];
    for (int i = 0; i < m; i++)
        delete[] M[i];
    delete[] M;
    return;
}

//---------------------------------------------------------------------------
//SET TO ZERO
//---------------------------------------------------------------------------

template <class Type>
Type** setToZeroMatrix(Type** M, int lx, int ly)
{
    for (int x=0 ; x<lx ; x++)
        for (int y=0 ; y<ly ; y++)
            M[x][y] = 0;
    return M;
}

template <class Type>
Type*** setToZeroMatrix(Type*** M, int lx, int ly, int k)
{
    for (int x=0 ; x<lx ; x++)
        for (int y=0 ; y<ly ; y++)
            for (int z=0 ; z<k ; z++)
                M[x][y][z] = 0;
    return M;
}

template <class Type>
Type** setToZeroMatrix(Type** M, int lx_min, int lx_max, int ly_min, int ly_max)
{
    for (int x=lx_min ; x<lx_max ; x++)
        for (int y=ly_min ; y<ly_max ; y++)
            M[x][y] = 0;
    return M;
}

template <class Type>
Type*** setToZeroMatrix(Type*** M, int lx_min, int lx_max, int ly_min, int ly_max, int k)
{
    for (int x=lx_min ; x<lx_max ; x++)
        for (int y=ly_min ; y<ly_max ; y++)
            for (int z=0 ; z<k ; z++)
                M[x][y][z] = 0;
    return M;
}

//---------------------------------------------------------------------------
//COPY
//---------------------------------------------------------------------------

template <class Type>
void copyMatrix(Type** M_Destinatario, Type** M_Sorgente, int lx, int ly)
{
    for (int x=0 ; x<lx ; x++)
        for (int y=0 ; y<ly ; y++)
            if ( M_Destinatario[x][y] != M_Sorgente[x][y] )
            M_Destinatario[x][y] = M_Sorgente[x][y];
}

template <class Type>
void copyMatrix(Type** M_Destinatario, Type** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max)
{
    for (int x=lx_min ; x<lx_max ; x++)
        for (int y=ly_min ; y<ly_max ; y++)
            if ( M_Destinatario[x][y] != M_Sorgente[x][y] )
            M_Destinatario[x][y] = M_Sorgente[x][y];
}

//---------------------------------------------------------------------------
//addition and subtraction
//---------------------------------------------------------------------------

template <class Type>
void addMatrices(Type **Dest, Type** M1, Type** M2, int lx, int ly)
{
    for (int x=0 ; x<lx ; x++)
        for (int y=0 ; y<ly ; y++)
            Dest[x][y] = M1[x][y] + M2[x][y];
    return;
}

template <class Type>
void subtractMatrix(Type **Dest, Type** M1, Type** menoM2, int lx, int ly)
{
    for (int x=0 ; x<lx ; x++)
        for (int y=0 ; y<ly ; y++)
            Dest[x][y] = M1[x][y] - menoM2[x][y];
    return;
}

//---------------------------------------------------------------------------
//INPUT - OUTPUT
//---------------------------------------------------------------------------

bool** readMatrix(bool** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
        for (int x=0 ; x<lx ; x++)
        {
            fscanf(f,"%s",&str);
			M[x][y] = atoi(str) == 0 ? false: true;
        }
    return M;
}

int** readMatrix(int** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
        for (int x=0 ; x<lx ; x++)
        {
            fscanf(f,"%s",&str);
            M[x][y] = atoi(str);
        }
    return M;
}

float** readMatrix(float** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
        for (int x=0 ; x<lx ; x++)
        {
            fscanf(f,"%s",&str);
            M[x][y] = (float)atof(str);
        }
    return M;
}

double** readMatrix(double** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
        for (int x=0 ; x<lx ; x++)
        {
            fscanf(f,"%s",&str);
            M[x][y] = atof(str);
        }
    return M;
}



void saveMatrix(bool** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            sprintf(str, "%d", M[x][y]);
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

void saveMatrix(int** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            sprintf(str, "%d", M[x][y]);
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

void saveMatrix(float** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            sprintf(str, "%.6f", M[x][y]);
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

void saveMatrix(double** M, int lx, int ly, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            sprintf(str, "%.6f", M[x][y]);
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

//*****************************************************

void saveMatrix(int** M, int lx, int ly, int threshold, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            if ( M[x][y] >= threshold )
                sprintf(str, "%d", M[x][y]);
            else
                strcpy(str, "0");
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

void saveMatrix(double** M, int lx, int ly, double threshold, FILE* f)
{
    char str[255];

    for (int y=0 ; y<ly ; y++)
    {
        for (int x=0 ; x<lx ; x++)
        {
            if ( M[x][y] >= threshold )
                sprintf(str, "%.6f", M[x][y]);
            else
                strcpy(str, "0.00");
            fprintf(f,"%s ",str);
        }
        fprintf(f,"\n");
    }
}

/***************************************************************************/
template <class Type> Type** readMatrix(Type** M, int lx, int ly, char* path)
{
	FILE *f;
	f = fopen(path, "r");
	M = readMatrix(M, lx, ly, f);
	fclose(f);
	return M;
}

template <class Type> void saveMatrix(Type** M, int lx, int ly, char* path)
{
	FILE *f;
	f = fopen(path, "w");
	saveMatrix(M, lx, ly, f);
	fclose(f);
}

template <class Type> void saveMatrix(Type** M, int lx, int ly, Type threshold, char* path)
{
	FILE *f;
	f = fopen(path, "w");
	saveMatrix(M, lx, ly, threshold, f);
	fclose(f);
}

/***************************************************************************/
//---------------------------------------------------------------------------
//DEFINITION OF TEMPLATE FUNCTIONS
//---------------------------------------------------------------------------

template bool** allocateMatrix(bool** M, int m, int n);
template int** allocateMatrix(int** M, int m, int n);
template float** allocateMatrix(float** M, int m, int n);
template double** allocateMatrix(double** M, int m, int n);

template double*** allocateMatrix(double*** M, int m, int n, int k);
template void deAllocateMatrix(double*** M, int m, int n, int k);

template void deAllocateMatrix(bool** M, int m, int n);
template void deAllocateMatrix(int** M, int m, int n);
template void deAllocateMatrix(float** M, int m, int n);
template void deAllocateMatrix(double** M, int m, int n);


template bool** setToZeroMatrix(bool** M, int lx, int ly);
template int** setToZeroMatrix(int** M, int lx, int ly);
template float** setToZeroMatrix(float** M, int lx, int ly);
template double** setToZeroMatrix(double** M, int lx, int ly);

template double*** setToZeroMatrix(double*** M, int lx, int ly, int k);

template bool** setToZeroMatrix(bool** M, int lx_min, int lx_max, int ly_min, int ly_max);
template int** setToZeroMatrix(int** M, int lx_min, int lx_max, int ly_min, int ly_max);
template float** setToZeroMatrix(float** M, int lx_min, int lx_max, int ly_min, int ly_max);
template double** setToZeroMatrix(double** M, int lx_min, int lx_max, int ly_min, int ly_max);

template double*** setToZeroMatrix(double*** M, int lx_min, int lx_max, int ly_min, int ly_max, int k);

template void copyMatrix(bool** M_Destinatario, bool** M_Sorgente, int lx, int ly);
template void copyMatrix(int** M_Destinatario, int** M_Sorgente, int lx, int ly);
template void copyMatrix(float** M_Destinatario, float** M_Sorgente, int lx, int ly);
template void copyMatrix(double** M_Destinatario, double** M_Sorgente, int lx, int ly);

template void copyMatrix(bool** M_Destinatario, bool** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max);
template void copyMatrix(int** M_Destinatario, int** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max);
template void copyMatrix(float** M_Destinatario, float** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max);
template void copyMatrix(double** M_Destinatario, double** M_Sorgente, int lx_min, int lx_max, int ly_min, int ly_max);

template void subtractMatrix(bool **Dest, bool** M1, bool** menoM2, int lx, int ly);
template void subtractMatrix(int **Dest, int** M1, int** menoM2, int lx, int ly);
template void subtractMatrix(float **Dest, float** M1, float** menoM2, int lx, int ly);
template void subtractMatrix(double **Dest, double** M1, double** menoM2, int lx, int ly);

template void addMatrices(bool **Dest, bool** M1, bool** M2, int lx, int ly);
template void addMatrices(int **Dest, int** M1, int** M2, int lx, int ly);
template void addMatrices(float **Dest, float** M1, float** M2, int lx, int ly);
template void addMatrices(double **Dest, double** M1, double** M2, int lx, int ly);

template bool** readMatrix(bool** M, int lx, int ly, char* path);
template int** readMatrix(int** M, int lx, int ly, char* path);
template float** readMatrix(float** M, int lx, int ly, char* path);
template double** readMatrix(double** M, int lx, int ly, char* path);

template void saveMatrix(bool** M, int lx, int ly, char* path);
template void saveMatrix(int** M, int lx, int ly, char* path);
template void saveMatrix(float** M, int lx, int ly, char* path);
template void saveMatrix(double** M, int lx, int ly, char* path);

template void saveMatrix(int** M, int lx, int ly, int threshold, char* path);
template void saveMatrix(double** M, int lx, int ly, double threshold, char* path);
