#include <stdio.h>
#include <math.h>

#define MAX_ITER 1000
#define EPSILON 0.0001
#define N 4  

// Функция для проверки условия сходимости (Диагональное преобладание)
int check_convergence(double A[N][N]) {
    int row_dominance = 1;  
    int col_dominance = 1;  
    
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                row_sum += fabs(A[i][j]);
            }
        }
        if (fabs(A[i][i]) <= row_sum) {
            row_dominance = 0;
            break;
        }
    }
    
    for (int j = 0; j < N; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < N; i++) {
            if (i != j) {
                col_sum += fabs(A[i][j]);
            }
        }
        if (fabs(A[j][j]) <= col_sum) {
            col_dominance = 0;
            break;
        }
    }
    
    return (row_dominance && col_dominance); // ||
}

int simple_iteration(double A[N][N], double b[N], double x[N]) {
    double x_old[N];
    double alpha[N]; //b
    double alpha_mat[N][N]; //a
    
    if (!check_convergence(A)) {
        printf("Условие сходимости не выполняется!\n");
        return -1;
    }
    
    // Преобразование системы к виду x = alpha*x + beta
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                alpha_mat[i][j] = 0.0;
            } else {
                alpha_mat[i][j] = -A[i][j] / A[i][i];
            }
        }
        alpha[i] = b[i] / A[i][i];
        x[i] = alpha[i]; // Начальное приближение
    }
    
    int iter = 0;
    double max_diff;
    
    printf("Iter\tx1\t\tx2\t\tx3\t\tx4\t\tMax Error\n");
    printf("----------------------------------------------------------------------------------------\n");
    
    do {
        // Сохраняем предыдущее приближение
        for (int i = 0; i < N; i++) {
            x_old[i] = x[i];
        }
        
        // Вычисляем новое приближение x⁽ᵏ⁺¹⁾ = αx⁽ᵏ⁾ + β
        for (int i = 0; i < N; i++) {
            x[i] = alpha[i];
            for (int j = 0; j < N; j++) {
                x[i] += alpha_mat[i][j] * x_old[j];
            }
        }
        
        max_diff = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = fabs(x[i] - x_old[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        
        iter++;
        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               iter, x[0], x[1], x[2], x[3], max_diff);
        
    } while (max_diff > EPSILON && iter < MAX_ITER);
    
    return iter;
}

int main() {
    double A[N][N] = {
        {15,  4,   -8,   2},
        {-3, 17,   -4,   7},
        {6,   5,  -16,   3},
        {3,  -7,    4,  -15}
    };
    
    double b[N] = {25, 14, 23, 30};
    double x[N];
    
    printf("Решение системы уравнений методом простой итерации\n");
    printf("Точность: %f\n\n", EPSILON);
    
    int iterations = simple_iteration(A, b, x);
    
    if (iterations > 0) {
        printf("\nРезультат:\n");
        printf("Количество итераций: %d\n", iterations);
        printf("x1 = %.6f\n", x[0]);
        printf("x2 = %.6f\n", x[1]);
        printf("x3 = %.6f\n", x[2]);
        printf("x4 = %.6f\n", x[3]);
        
    }
    
    return 0;
}