#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int UL_decomposition(int n, double A[n][n], double LU[n][n]) {
    // копируем исходную матрицу в LU
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            LU[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
        // Вычисляем элементы L для столбца i (ниже диагонали)
        for (int j = i + 1; j < n; j++) {
            if (fabs(LU[i][i]) < 1e-12) {
                printf("Ошибка: нулевой диагональный элемент U[%d][%d]\n", i, i);
                return 0;
            }
            LU[j][i] /= LU[i][i]; 
            
            for (int k = i + 1; k < n; k++) {
                LU[j][k] -= LU[j][i] * LU[i][k];
            }
        }
    }
    return 1;
}

// Ly = b
void solve_L(int n, double LU[n][n], double b[n], double y[n]) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];  // LU[i][j] содержит L[i][j]
        }
        // L[i][i] = 1, поэтому не делим
    }
}

// Ux = y (с верхней треугольной матрицей U)
void solve_U(int n, double LU[n][n], double y[n], double x[n]) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];  // LU[i][j] содержит U[i][j]
        }
        x[i] /= LU[i][i];  // LU[i][i] содержит U[i][i]
    }
}

void extract_L(int n, double LU[n][n], double L[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i > j) {
                L[i][j] = LU[i][j]; 
            } else if (i == j) {
                L[i][j] = 1.0;    
            } else {
                L[i][j] = 0.0;     
            }
        }
    }
}

void extract_U(int n, double LU[n][n], double U[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                U[i][j] = LU[i][j];  
            } else {
                U[i][j] = 0.0;      
            }
        }
    }
}

double determinant(int n, double LU[n][n]) {
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= LU[i][i]; 
    }
    return det;
}

int inverse_matrix(int n, double A[n][n], double invA[n][n]) {
    double LU[n][n];
    double y[n], x[n];
    double e[n];
    
    if (!UL_decomposition(n, A, LU)) {
        return 0;
    }
    
    // для каждого столбца единичной матрицы решаем систему
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            e[i] = (i == j) ? 1.0 : 0.0;
        }
        
        // Ly = e
        solve_L(n, LU, e, y);
        // Ux = y
        solve_U(n, LU, y, x);
        
        // записываем результат в j-ый столбец обратной матрицы
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

void print_LU(int n, double LU[n][n]) {
    printf("Объединенная матрица LU:\n");
    printf("(L-элементы ниже диагонали, U-элементы на и выше диагонали)\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i > j) {
                printf("L%10.6f ", LU[i][j]);
            } else {
                printf("U%10.6f ", LU[i][j]);
            }
        }
        printf("\n");
    }
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
    
    printf("Исходная матрица:\n");
    print_matrix(n, n, A);
    
    printf("\nВектор правой части b:\n");
    print_vector(n, b);
    
    // Используем объединенную матрицу LU
    double LU[n][n];
    if (!UL_decomposition(n, A, LU)) {
        printf("Ошибка UL-разложения\n");
        return 1;
    }
    
    print_LU(n, LU);
    
    // можно извлечь L и U отдельно
    double L[n][n], U[n][n];
    extract_L(n, LU, L);
    extract_U(n, LU, U);
    
    printf("\nМатрица L (извлеченная):\n");
    print_matrix(n, n, L);
    
    printf("\nМатрица U (извлеченная):\n");
    print_matrix(n, n, U);
    
    // проверка: L * U = A
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
    
    // решение
    double y[n], x[n];
    solve_L(n, LU, b, y);
    solve_U(n, LU, y, x);
    
    printf("\nРешение системы Ax = b:\n");
    print_vector(n, x);
    
    double det = determinant(n, LU);
    printf("\nОпределитель: %.6f\n", det);
    
    double invA[n][n];
    if (inverse_matrix(n, A, invA)) {
        printf("\nОбратная матрица:\n");
        print_matrix(n, n, invA);
    }
    
    // Матрица 5x5
    printf("\n=== Матрица 5x5 ===\n");
    
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
    
    double LU5[n5][n5];
    if (UL_decomposition(n5, A5, LU5)) {
        double y5[n5], x5[n5];
        solve_L(n5, LU5, b5, y5);
        solve_U(n5, LU5, y5, x5);
        
        printf("\nРешение для матрицы 5x5:\n");
        print_vector(n5, x5);
        
        double det5 = determinant(n5, LU5);
        printf("\nОпределитель: %.6f\n", det5);

        double invA5[n5][n5];
        if (inverse_matrix(n5, A5, invA5)) {
            printf("\nОбратная матрица для 5x5:\n");
            print_matrix(n5, n5, invA5);
        }
    }
    
    return 0;
}