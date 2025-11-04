#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Функция для вычисления ПЕРВОЙ ПРОИЗВОДНОЙ сплайна в точке xp
double spline_derivative(int n, double x[], double a[], double b[], double c[], double d[], double xp) {
    for (int i = 0; i < n; i++) {
        if (xp >= x[i] && xp <= x[i+1]) {
            double dx = xp - x[i];
            // Производная полинома: S'(x) = b_i + 2*c_i*dx + 3*d_i*dx^2
            return b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;
        }
    }
    return 0.0;
}

//метод прогонки
void classical_sweep_method(int n, double a[], double b[], double c[], double d_vec[], double x[]) {
    double *P = (double*)malloc(n * sizeof(double));
    double *Q = (double*)malloc(n * sizeof(double));
    
    P[0] = -c[0] / b[0];
    Q[0] = d_vec[0] / b[0];
    
    for (int i = 1; i < n-1; i++) {
        double denominator = b[i] + a[i] * P[i-1];
        P[i] = -c[i] / denominator;
        Q[i] = (d_vec[i] - a[i] * Q[i-1]) / denominator;
    }
    
    int last = n-1;
    double denominator = b[last] + a[last] * P[last-1];
    P[last] = 0.0;  
    Q[last] = (d_vec[last] - a[last] * Q[last-1]) / denominator;
    
    x[last] = Q[last];
    for (int i = n-2; i >= 0; i--) {
        x[i] = P[i] * x[i+1] + Q[i];
    }
    
    free(P);
    free(Q);
}

void natural_cubic_spline(int n, double x[], double y[], double a[], double b[], double c[], double d[]) {
    double *h = (double*)malloc(n * sizeof(double));
    
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    c[0] = 0.0;
    c[n] = 0.0;    
    
    int system_size = n - 1; //всего n-1 неизвестных, тк c[0] = 0.0 и c[n] = 0.0
    
    double *A = (double*)malloc(system_size * sizeof(double)); // нижняя диагональ
    double *B = (double*)malloc(system_size * sizeof(double)); // главная диагональ
    double *C = (double*)malloc(system_size * sizeof(double)); // верхняя диагональ  
    double *D = (double*)malloc(system_size * sizeof(double)); // правая часть
    double *c_solution = (double*)malloc(system_size * sizeof(double)); // решения
    
    //hᵢ₋₁·cᵢ₋₁ + 2(hᵢ₋₁ + hᵢ)·cᵢ + hᵢ·cᵢ₊₁ = 3[(yᵢ₊₁ - yᵢ)/hᵢ - (yᵢ - yᵢ₋₁)/hᵢ₋₁]
    for (int i = 1; i < n; i++) {
        int sys_idx = i - 1; 
        
        A[sys_idx] = h[i-1];                    
        B[sys_idx] = 2.0 * (h[i-1] + h[i]);     
        C[sys_idx] = h[i];                    
        D[sys_idx] = 3.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
    }
    
    A[0] = 0.0;         
    C[system_size-1] = 0.0; 
    
    classical_sweep_method(system_size, A, B, C, D, c_solution);
    
    for (int i = 0; i < system_size; i++) {
        c[i + 1] = c_solution[i];
    }
    
    for (int i = 0; i < n; i++) {
        a[i] = y[i];
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
    }
    
    free(h);
    free(A);
    free(B);
    free(C);
    free(D);
    free(c_solution);
}

double spline_value(int n, double x[], double a[], double b[], double c[], double d[], double xp) {
    for (int i = 0; i < n; i++) {
        if (xp >= x[i] && xp <= x[i+1]) {
            double dx = xp - x[i];
            return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    }
    return 0.0; 
}

void check_continuity(int n, double x[], double a[], double b[], double c[], double d[]) {
    printf("\n=== ПРОВЕРКА НЕПРЕРЫВНОСТИ ПРОИЗВОДНЫХ ===\n");
    
    for (int i = 1; i < n; i++) {
        // Левая производная в узле x[i] (от отрезка i-1)
        double dx_left = x[i] - x[i-1];
        double S_left = b[i-1] + 2.0 * c[i-1] * dx_left + 3.0 * d[i-1] * dx_left * dx_left;
        
        // Правая производная в узле x[i] (от отрезка i)
        double S_right = b[i];
        
        // Вторая производная слева
        double S2_left = 2.0 * c[i-1] + 6.0 * d[i-1] * dx_left;
        
        // Вторая производная справа
        double S2_right = 2.0 * c[i];
        
        printf("Узел x[%d] = %.3f:\n", i, x[i]);
        printf("  S'(слева) = %.6f, S'(справа) = %.6f, разность = %.6f\n", 
               S_left, S_right, fabs(S_left - S_right));
        printf("  S''(слева) = %.6f, S''(справа) = %.6f, разность = %.6f\n", 
               S2_left, S2_right, fabs(S2_left - S2_right));
    }
}

int main() {
    int n = 10; //количество отрезков
    double x[] = {0.06, 0.356, 0.763, 1.096, 1.54, 2.021, 2.354, 2.761, 3.205, 3.464, 3.76};
    double y[] = {2.478, 3.471, 3.257, 1.086, 1.291, 3.735, 0.729, 0.684, 2.371, 3.142, 2.846};
    
    double *a = (double*)malloc(n * sizeof(double));
    double *b = (double*)malloc(n * sizeof(double));
    double *c = (double*)malloc((n+1) * sizeof(double));
    double *d = (double*)malloc(n * sizeof(double));
    
    natural_cubic_spline(n, x, y, a, b, c, d);
    
    double x_star = 1.285;
    
    //находим отрезок для x_star
    int segment = -1;
    for (int i = 0; i < n; i++) {
        if (x_star >= x[i] && x_star <= x[i+1]) {
            segment = i;
            break;
        }
    }
    
    if (segment == -1) {
        printf("Точка x* = %.3f вне диапазона данных\n", x_star);
        return 1;
    }
    
    //вычисляем значение в точке x_star
    double y_star = spline_value(n, x, a, b, c, d, x_star);
    
    printf("=== ЕСТЕСТВЕННЫЙ КУБИЧЕСКИЙ СПЛАЙН ДЕФЕКТА 1 ===\n");
    printf("Точка x* = %.3f находится на отрезке [%.3f, %.3f] (i = %d)\n", 
           x_star, x[segment], x[segment+1], segment+1);
    printf("Коэффициенты сплайна на этом отрезке:\n");
    printf("a[%d] = %.6f\n", segment+1, a[segment]);
    printf("b[%d] = %.6f\n", segment+1, b[segment]);
    printf("c[%d] = %.6f\n", segment+1, c[segment]);
    printf("d[%d] = %.6f\n", segment+1, d[segment]);
    printf("\n");
    printf("Значение сплайна в точке x* = %.3f:\n", x_star);
    printf("S(%.3f) = %.6f\n", x_star, y_star);
    printf("\n");
    
    printf("ВСЕ КОЭФФИЦИЕНТЫ СПЛАЙНА:\n");
    printf("i\tx[i]\ty[i]\ta[i]\t\tb[i]\t\tc[i]\t\td[i]\n");
    for (int i = 0; i < n; i++) {
        printf("%d\t%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               i, x[i], y[i], a[i], b[i], c[i], d[i]);
    }
    printf("%d\t%.3f\t%.3f\t\t\t\t\t%.6f\n", n, x[n], y[n], c[n]);

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}