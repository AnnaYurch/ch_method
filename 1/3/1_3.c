#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8

int main() {
    
    /*
    double A[N][N] = {
        {3, -1,  0,  0,  0,  0,  0,  0},
        {3,  7, -2,  0,  0,  0,  0,  0},
        {0,  2,  3, -1,  0,  0,  0,  0},
        {0,  0,  1, -4, -3,  0,  0,  0},
        {0,  0,  0,  3,  8,  1,  0,  0},
        {0,  0,  0,  0, -5, 10,  4,  0},
        {0,  0,  0,  0,  0,  3,  5, -2},
        {0,  0,  0,  0,  0,  0,  2,  3}
    };

    double d[N] = {15, -13, 1, -3, 22, 3, 7, 13};
    */

    double A[N][N] = {
        {12, 5,  0,  0,  0,  0,  0,  0},
        {5,  10, -4, 0,  0,  0,  0,  0},
        {0,  4,  9,  -5, 0,  0,  0,  0},
        {0,  0,  7,  11, 4,  0,  0,  0},
        {0,  0,  0,  6,  10, 3,  0,  0},
        {0,  0,  0,  0,  -5, 13, 7,  0},
        {0,  0,  0,  0,  0,  5,  9,  -2},
        {0,  0,  0,  0,  0,  0,  -8, 9}
        };

    double d[N+1] = {0, -71, -69, -54, -54, -19, -17, -48, 76};
    
    double a[N+1], b[N+1], c[N+1];
    double P[N+1], Q[N+1];
    double x[N+1];
    
    for (int i = 1; i <= N; i++) {
        b[i] = A[i-1][i-1];
        
        if (i > 1) {
            a[i] = A[i-1][i-2];
        } else {
            a[i] = 0;
        }
        
        if (i < N) {
            c[i] = A[i-1][i];
        } else {
            c[i] = 0;
        }
    }

    printf("Коэффициенты a, b, c, d:\n");
    printf(" i | a[i] | b[i] | c[i] | d[i] \n");
    printf("---+------+------+------+------\n");
    for (int i = 1; i <= N; i++) {
        printf("%2d | %4.0f | %4.0f | %4.0f | %4.0f\n", 
               i, a[i], b[i], c[i], d[i-1]);
    }

    printf("\nУстойчивость метода:\n");
    int stable = 1;
    for (int i = 1; i <= N; i++) {
        double left_side = fabs(b[i]);
        double right_side = fabs(a[i]) + fabs(c[i]);
        printf("|b[%d]| = %.1f, |a[%d]| + |c[%d]| = %.1f", 
               i, left_side, i, i, right_side);
        if (left_side >= right_side) {
            printf(" - выполняется\n");
        } else {
            printf(" - НЕ выполняется\n");
            stable = 0;
        } 
    }

    //краевые
    double c1_b1 = fabs(c[1] / b[1]);
    double an_bn = fabs(a[N] / b[N]);
    int edge1_ok = (c1_b1 < 1);
    int edge2_ok = (an_bn < 1);
    
    printf("\nКраевые условия:\n");
    printf("|c1/b1| = %.3f < 1: %s\n", c1_b1, edge1_ok ? "выполняется" : "НЕ выполняется");
    printf("|an/bn| = %.3f < 1: %s\n", an_bn, edge2_ok ? "выполняется" : "НЕ выполняется");
    
    if (stable && edge1_ok && edge2_ok) {
        printf("\nЗначит устойчив\n");
    } else {
        printf("\nЗначит неустойчив, лучше использовать другой метод\n");
        return 1;
    }
    
    printf("\nПрямой ход:\n");
    
    P[1] = -c[1] / b[1];
    Q[1] = d[1] / b[1];
    printf("P[1] = %.6f, Q[1] = %.6f\n", P[1], Q[1]);
    
    for (int i = 2; i <= N; i++) {
        double denominator = b[i] + a[i] * P[i-1];
        P[i] = -c[i] / denominator;
        Q[i] = (d[i] - a[i] * Q[i-1]) / denominator;
        printf("P[%d] = %.6f, Q[%d] = %.6f\n", i, P[i], i, Q[i]);
    }
    
    printf("\nОбратный ход:\n");
    
    x[N] = Q[N];
    printf("x[%d] = %.6f\n", N, x[N]);
    
    for (int i = N-1; i >= 1; i--) {
        x[i] = P[i] * x[i+1] + Q[i];
        printf("x[%d] = %.6f\n", i, x[i]);
    }
    
    double determinant = 1.0;
    for (int i = 1; i <= N; i++) {
        if (i == 1) {
            determinant *= b[1];
        } else {
            determinant *= (b[i] + a[i] * P[i-1]);
        }
    }
    
    printf("\nОпределитель: %.0f\n", determinant);
    
    return 0;
}