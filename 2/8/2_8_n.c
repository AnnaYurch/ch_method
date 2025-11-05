#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS 0.0001
#define MAX_ITER 1000

int UL_decomposition(int n, double A[n][n], double LU[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            LU[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
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
            y[i] -= LU[i][j] * y[j];
        }
    }
}

// Ux = y
void solve_U(int n, double LU[n][n], double y[n], double x[n]) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= LU[i][i];
    }
}

int inverse_matrix(int n, double A[n][n], double invA[n][n]) {
    double LU[n][n];
    double y[n], x[n];
    double e[n];
    
    if (!UL_decomposition(n, A, LU)) {
        return 0;
    }
    
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            e[i] = (i == j) ? 1.0 : 0.0;
        }
        
        solve_L(n, LU, e, y);
        solve_U(n, LU, y, x);
        
        for (int i = 0; i < n; i++) {
            invA[i][j] = x[i];
        }
    }
    return 1;
}

void f(double x1, double x2, double *f1, double *f2) {
    *f1 = x1*x1 + x2*x2 - 5*sin(x1) - 9;
    *f2 = 2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3;
}

void jacobian(double x1, double x2, double J[2][2]) {
    J[0][0] = 2*x1 - 5*cos(x1);
    J[0][1] = 2*x2;
    J[1][0] = 4*x1 + 2*x2 - 4 + x2*sin(x1);
    J[1][1] = 2*x1 - 6*x2 - cos(x1);
}

int newton_method(double *x1, double *x2) {
    double x1_old, x2_old;
    double f1, f2;
    double J[2][2], J_inv[2][2];
    double dx[2];
    int iter = 0;
    double dx_norm = 1.0;
    
    printf("Метод Ньютона:\n");
    printf("k\tx1\t\tx2\t\tf1\t\tf2\n");
    
    do {
        x1_old = *x1;
        x2_old = *x2;
        
        f(*x1, *x2, &f1, &f2);
        
        if (fabs(f1) < EPS && fabs(f2) < EPS) break;
        
        jacobian(*x1, *x2, J);
        
        if (!inverse_matrix(2, J, J_inv)) {
            printf("ОШИБКА: Определитель матрицы Якоби равен нулю!\n");
            return -1; 
        }
        
        // ΔX = -J⁻¹·F
        dx[0] = -(J_inv[0][0] * f1 + J_inv[0][1] * f2);
        dx[1] = -(J_inv[1][0] * f1 + J_inv[1][1] * f2);
        
        *x1 += dx[0];
        *x2 += dx[1];
        
        dx_norm = sqrt((*x1 - x1_old)*(*x1 - x1_old) + 
                      (*x2 - x2_old)*(*x2 - x2_old));
        
        printf("%d\t%.6f\t%.6f\t%.2e\t%.2e\n", iter, *x1, *x2, f1, f2);
        
        if (++iter >= MAX_ITER) break;
        
    } while (dx_norm > EPS);
    
    return iter;
}

//==============

int main() {
    /*
        x1=-1.318459, x2=1.555637 
        x1=-0.881570, x2=-2.089052 
        x1=2.423889, x2=2.532399 
        x1=2.987711, x2=-0.916494 
    */
    
    double solutions[][2] = {
        {2.0, 2.5},   
        {-1.0, -2.0},  
        {-1.5, 1.5},   
        {1.3, -1.8}
    };
    
    for (int i = 0; i < 4; i++) { 
        printf("\n=== Корень %d ===\n", i+1);
        
        double x1_n = solutions[i][0], x2_n = solutions[i][1];
        
        printf("\n1. ");
        int iter_n = newton_method(&x1_n, &x2_n);
        
        printf("\nИтоговые решения:\n");
        printf("Ньютон:    x1=%.6f, x2=%.6f (итераций: %d)\n", x1_n, x2_n, iter_n);
    }
    
    return 0;
}