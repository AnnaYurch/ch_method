#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5
#define EXTENDED_N (2 * N + 1)

void printExtendedMatrix(double mat[N][EXTENDED_N]) {
    printf("Расширенная матрица [A|B|A^(-1)]:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%8.4f ", mat[i][j]);
        }
        printf(" | ");
        
        printf("%8.4f | ", mat[i][N]);
        
        for (int j = N + 1; j < EXTENDED_N; j++) {
            printf("%8.4f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void gaussJordanExtended(double mat[N][EXTENDED_N], double x[N], double inv[N][N]) {
    double temp;
    int i, j, k;
    
    for (k = 0; k < N; k++) {
        //нормализуем текущую строку
        temp = mat[k][k];
        for (j = k; j < EXTENDED_N; j++) {
            mat[k][j] /= temp;
        }
        
        //исключаем переменную из других строк
        for (i = 0; i < N; i++) {
            if (i != k) {
                temp = mat[i][k];
                for (j = k; j < EXTENDED_N; j++) {
                    mat[i][j] -= temp * mat[k][j]; 
                }
            }
        }
        printExtendedMatrix(mat);
    }
    
    //извлекаем решение и обратную матрицу
    for (i = 0; i < N; i++) {
        x[i] = mat[i][N];
        
        for (j = 0; j < N; j++) {
            inv[i][j] = mat[i][N + 1 + j];
        }
    }
}

double determinant(double A[N][N]) {
    double det = 1.0;
    double temp;
    double B[N][N];
    int i, j, k;
    
    //копируем матрицу
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            B[i][j] = A[i][j];
        }
    }
    
    //приводим к треугольному виду
    for (k = 0; k < N-1; k++) {
        for (i = k+1; i < N; i++) {
            temp = B[i][k] / B[k][k];
            for (j = k; j < N; j++) {
                B[i][j] -= temp * B[k][j];
            }
        }
    }
    
    //вычисляем определитель
    for (i = 0; i < N; i++) {
        det *= B[i][i];
    }
    
    return det;
}

int main() {
    double extendedMat[N][EXTENDED_N] = {0};
    
    /*
    double A[N][N] = {
        {2, -3, -4, -1},
        {3, 2, 1, -5},
        {-6, 3, -5, 2},
        {7, 5, -4, 2}
    };
    
    double b[N] = {15, -12, 10, 18};
    */
   double A[N][N] = {
        {2, -3, -4, 3, 7},
        {4, 3, 5, -4, 2},
        {-2, 3, -5, 2, -3},
        {6, -1, -4, 1, 3},
        {-3, 9, -5, 4, -2}
    };

    double b[N] = {-17, -20, -11, -17, -14};
    //делаем расширенную матрицу
    for (int i = 0; i < N; i++) {
        //копируем матрицу A
        for (int j = 0; j < N; j++) {
            extendedMat[i][j] = A[i][j];
        }
        
        //копируем вектор B
        extendedMat[i][N] = b[i];
        
        //добавляем единичную матрицу
        for (int j = N + 1; j < EXTENDED_N; j++) {
            extendedMat[i][j] = (j == N + 1 + i) ? 1.0 : 0.0;
        }
    }
    
    printf("Исходная расширенная матрица:\n");
    printExtendedMatrix(extendedMat);
    
    double x[N];  //решение
    double inv[N][N];  //обратная матрица
    
    gaussJordanExtended(extendedMat, x, inv);
    
    printf("Решение системы:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %8.4f\n", i+1, x[i]);
    }
    printf("\n");
    
    double det = determinant(A);
    printf("Определитель матрицы: %.4f\n\n", det);
    
    printf("Обратная матрица:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%8.4f ", inv[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}