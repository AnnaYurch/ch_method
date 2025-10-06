#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define N 4
#define EPS 0.0001
#define MAX_ITER 100

bool checkConvergence(double A[N][N]) {
    bool row_dominance = true;
    bool col_dominance = true;
    
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) sum += fabs(A[i][j]);
        }
        if (fabs(A[i][i]) <= sum) {
            row_dominance = false;
        }
    }
    
    for (int j = 0; j < N; j++) {
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            if (i != j) sum += fabs(A[i][j]);
        }
        if (fabs(A[j][j]) <= sum) {
            col_dominance = false;
        }
    }
    
    return (row_dominance || col_dominance); 

}

//метод Зейделя
void seidelMethod(double A[N][N], double b[N], double x[N], int *iterations) {
    double x_old[N];
    double error;
    int iter = 0;
    
    printf("\nIteration\tx1\t\tx2\t\tx3\t\tx4\t\tMax Error\n");
    
    do {
        for (int i = 0; i < N; i++) {
            x_old[i] = x[i];
        }
        
        for (int i = 0; i < N; i++) {
            double sum = b[i];
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    sum -= A[i][j] * x[j];
                }
            }
            x[i] = sum / A[i][i];
        }
        
        error = 0.0;
        for (int i = 0; i < N; i++) {
            double current_error = fabs(x[i] - x_old[i]);
            if (current_error > error) {
                error = current_error;
            }
        }
        
        iter++;
        printf("%d\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               iter, x[0], x[1], x[2], x[3], error);
        
    } while (error > EPS && iter < MAX_ITER);
    
    *iterations = iter;
}

int main() {
    
    double A[N][N] = {
        {15,  4,   -8,   2},
        {-3, 17,   -4,   7},
        {6,   5,  -16,   3},
        {3,  -7,    4,  -15}
    };
    
    double b[N] = {25, 14, 23, 30};
    
    double x[N] = {0, 0, 0, 0};
    
    printf("Матрица A:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%6.1f ", A[i][j]);
        }
        printf(" | %6.1f\n", b[i]);
    }
    
    if (!checkConvergence(A)) {
        printf("Условие сходимости не выполняется!\n");
        return 1;
    }
    
    int iterations;
    seidelMethod(A, b, x, &iterations);
    
    printf("\nРезультат:\n");
    printf("Количество итераций: %d\n", iterations);
    printf("x1 = %.6f\n", x[0]);
    printf("x2 = %.6f\n", x[1]);
    printf("x3 = %.6f\n", x[2]);
    printf("x4 = %.6f\n", x[3]);
    
    return 0;
}