#include <stdio.h>
#include <stdlib.h>

void natural_cubic_spline(int n, double x[], double y[], double a[], double b[], double c[], double d[]) {
    double *h = (double*)malloc(n * sizeof(double)); // длина i-го отрезка между точками
    double *alpha = (double*)malloc(n * sizeof(double));
    double *l = (double*)malloc((n+1) * sizeof(double));
    double *mu = (double*)malloc((n+1) * sizeof(double));
    double *z = (double*)malloc((n+1) * sizeof(double));
    
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    c[0] = 0.0;
    c[n] = 0.0;
    
    //решаем трёхдиагональную систему для c[i]
    for (int i = 1; i < n; i++) {
        alpha[i] = 3.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
    }
    
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    
    //метод прогонки
    for (int i = 1; i < n; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }
    
    l[n] = 1.0;
    z[n] = 0.0;
    
    for (int j = n-1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j+1];
    }
    
    for (int i = 0; i < n; i++) {
        a[i] = y[i];
        b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (c[i+1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i]);
    }
    
    free(h);
    free(alpha);
    free(l);
    free(mu);
    free(z);
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
    printf("i\tx[i]\ty[i]\ta[i]\tb[i]\tc[i]\td[i]\n");
    for (int i = 0; i < n; i++) {
        printf("%d\t%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               i, x[i], y[i], a[i], b[i], c[i], d[i]);
    }
    printf("%d\t%.3f\t%.3f\t\t\t%.6f\n", n, x[n], y[n], c[n]);
    
    free(a);
    free(b);
    free(c);
    free(d);
    
    return 0;
}