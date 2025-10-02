#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int UL_decomposition(int n, double A[n][n], double U[n][n], double L[n][n]) {
    //инициализация U и L
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            U[i][j] = 0.0;
            L[i][j] = 0.0;
        }
        L[i][i] = 1.0; //диагональ L = 1
    }
    
    //заполнение U и L
    for (int i = 0; i < n; i++) {
        //вычисляем элементы U для строки i
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
        
        // Вычисляем элементы L для столбца i
        for (int j = i + 1; j < n; j++) {
            if (fabs(U[i][i]) < 1e-12) {
                printf("Ошибка: нулевой диагональный элемент");
                return 0;
            }
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
    return 1;
}

//Ly = b
void solve_L(int n, double L[n][n], double b[n], double y[n]) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
}

//Ux = y
void solve_U(int n, double U[n][n], double y[n], double x[n]) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

double determinant(int n, double U[n][n]) {
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= U[i][i];
    }
    return det;
}

int inverse_matrix(int n, double A[n][n], double invA[n][n]) {
    double U[n][n], L[n][n];
    double y[n], x[n];
    double e[n];
    
    if (!UL_decomposition(n, A, U, L)) {
        return 0;
    }
    
    // Для каждого столбца единичной матрицы решаем систему
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            e[i] = (i == j) ? 1.0 : 0.0;
        }
        
        //Ly = e
        solve_L(n, L, e, y);
        //Ux = y
        solve_U(n, U, y, x);
        
        // Записываем результат в j-ый столбец обратной матрицы
        for (int i = 0; i < n; i++) {
            invA[i][j] = x[i];
        }
    }
    return 1;
}

void print_matrix(int n, int m, double matrix[n][m]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%10.6f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_vector(int n, double vector[n]) {
    for (int i = 0; i < n; i++) {
        printf("%10.6f ", vector[i]);
    }
    printf("\n");
}

int main() {
    int n = 4;
    
    double A[4][4] = {
        {2, -3, -4, -1},
        {3, 2, 1, -5},
        {-6, 3, -5, 2},
        {7, 5, -4, 2}
    };
    
    double b[4] = {15, -12, 10, 18};
    
    printf("Исходная матрица:");
    print_matrix(n, n, A);
    
    printf("\nВектор правой части b:\n");
    print_vector(n, b);
    
    double U[n][n], L[n][n];
    if (!UL_decomposition(n, A, U, L)) {
        printf("Ошибка UL-разложения\n");
        return 1;
    }
    
    printf("\nМатрица L:\n");
    print_matrix(n, n, L);
    
    printf("\nМатрица U:\n");
    print_matrix(n, n, U);
    
    //A = L * U
    double check[n][n];
    printf("\nПроверка: L * U =\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            check[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                check[i][j] += L[i][k] * U[k][j];
            }
        }
    }
    print_matrix(n, n, check);
    
    //решение системы
    double y[n], x[n];
    solve_L(n, L, b, y);
    solve_U(n, U, y, x);
    
    printf("\nРешение системы Ax = b:\n");
    print_vector(n, x);
    
    double det = determinant(n, U);
    printf("\nОпределитель: %.6f\n", det);
    
    double invA[n][n];
    if (inverse_matrix(n, A, invA)) {
        printf("\nОбратная матрица:\n");
        print_matrix(n, n, invA);
    }
    
    //моя матрица 5x5
    printf("\n=== матрица 5x5 ===\n");
    
    int n5 = 5;
    double A5[5][5] = {
        {2, -3, -4,  3,  7},
        {4,  3,  5, -4,  2},
        {-2, 3, -5,  2, -3},
        {6, -1, -4,  1,  3},
        {-3, 9, -5,  4, -2}
    };
    
    double b5[5] = {-17, -20, -11, -17, -14};
    
    printf("Матрица 5x5:\n");
    print_matrix(n5, n5, A5);
    
    double U5[n5][n5], L5[n5][n5];
    if (UL_decomposition(n5, A5, U5, L5)) {
        double y5[n5], x5[n5];
        solve_L(n5, L5, b5, y5);
        solve_U(n5, U5, y5, x5);
        
        printf("\nРешение для матрицы 5x5:\n");
        print_vector(n5, x5);
        
        double det5 = determinant(n5, U5);
        printf("\nОпределитель: %.6f\n", det5);

        double invA5[n5][n5];
        if (inverse_matrix(n5, A5, invA5)) {
            printf("\nОбратная матрица для 5x5:\n");
            print_matrix(n5, n5, invA5);
        }
    }
    
    return 0;
}