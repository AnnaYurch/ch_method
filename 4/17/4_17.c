#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NODES 100

// Функция метода Гаусса для решения СЛАУ
void gauss_solve(int n, double A[MAX_NODES][MAX_NODES], double b[MAX_NODES], double x[MAX_NODES]) {
    double augmented[MAX_NODES][MAX_NODES + 1];
    int i, j, k;
    
    // Создаем расширенную матрицу [A | b]
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = b[i];
    }
    
    // Прямой ход метода Гаусса
    for (k = 0; k < n; k++) {
        // Поиск главного элемента
        int max_row = k;
        for (i = k + 1; i < n; i++) {
            if (fabs(augmented[i][k]) > fabs(augmented[max_row][k])) {
                max_row = i;
            }
        }
        
        // Перестановка строк
        if (max_row != k) {
            for (j = k; j <= n; j++) {
                double temp = augmented[k][j];
                augmented[k][j] = augmented[max_row][j];
                augmented[max_row][j] = temp;
            }
        }
        
        // Нормализация текущей строки
        double pivot = augmented[k][k];
        for (j = k; j <= n; j++) {
            augmented[k][j] /= pivot;
        }
        
        // Исключение переменной из других строк
        for (i = 0; i < n; i++) {
            if (i != k) {
                double factor = augmented[i][k];
                for (j = k; j <= n; j++) {
                    augmented[i][j] -= factor * augmented[k][j];
                }
            }
        }
    }
    
    // Извлекаем решение
    for (i = 0; i < n; i++) {
        x[i] = augmented[i][n];
    }
}

double analytical_solution(double x) {
    return (2.0 * sin(x) + cos(x)) / (x * x - 1.0);
}

//коэффициенты 
double p(double x) {
    return (4.0 * x) / (x * x - 1.0);
}

double q(double x) {
    return (x * x + 1.0) / (x * x - 1.0);
}

double f(double x) {
    return 0.0;
}

void solve_bvp(int order, double h, int n, double x[MAX_NODES], double y[MAX_NODES]) {
    double A[MAX_NODES][MAX_NODES] = {0};
    double b[MAX_NODES] = {0};
    
    if (order == 1) {
        //1-го порядка
        
        // Левая граница: y'(-0.9) = -30.4
        A[0][0] = -1.0 / h; // коэффициент при y0
        A[0][1] = 1.0 / h; // коэффициент при y1
        b[0] = -30.4;
        
        // Правая граница: y'(0.9) = -57.9
        A[n-1][n-2] = -1.0 / h;
        A[n-1][n-1] = 1.0 / h;
        b[n-1] = -57.9;
        
    } else {
        //2-го порядка
        
        // Левая граница: (-3y₀ + 4y₁ - y₂)/(2h) = -30.4
        A[0][0] = -3.0 / (2.0 * h);
        A[0][1] = 4.0 / (2.0 * h);
        A[0][2] = -1.0 / (2.0 * h);
        b[0] = -30.4;
        
        // Правая граница: (y_{n-3} - 4y_{n-2} + 3y_{n-1})/(2h) = -57.9
        A[n-1][n-3] = 1.0 / (2.0 * h);
        A[n-1][n-2] = -4.0 / (2.0 * h);
        A[n-1][n-1] = 3.0 / (2.0 * h);
        b[n-1] = -57.9;
    }
    
    //внутренние точки
    for (int i = 1; i < n-1; i++) {
        double x_i = x[i];
        double p_i = p(x_i);
        double q_i = q(x_i);
        double f_i = f(x_i);
        
        // (y_{i+1} - 2y_i + y_{i-1})/h² + p_i*(y_{i+1} - y_{i-1})/(2h) + q_i*y_i = f_i
        
        A[i][i-1] = (1.0/(h*h)) - (p_i/(2.0*h));
        A[i][i]   = (-2.0/(h*h)) + q_i;
        A[i][i+1] = (1.0/(h*h)) + (p_i/(2.0*h));
        b[i] = f_i;
    }
    
    gauss_solve(n, A, b, y);
}

int main() {
    double a = -0.9, b = 0.9, h = 0.15;
    int n = (int)((b - a) / h) + 1;
    
    double x[MAX_NODES];
    double y1[MAX_NODES]; //решение с 1-м порядком
    double y2[MAX_NODES]; //решение со 2-м 
    double y_analytical[MAX_NODES];
    
    //создаем сетку
    for (int i = 0; i < n; i++) {
        x[i] = a + i * h;
    }
    
    printf("РЕШЕНИЕ КРАЕВОЙ ЗАДАЧИ №44\n");
    printf("Уравнение: (x² - 1)y'' + 4xy' + (x² + 1)y = 0\n");
    printf("Граничные условия: y'(-0.9) = -30.4, y'(0.9) = -57.9\n");
    printf("Шаг: h = %.2f, Количество узлов: %d\n\n", h, n);
    
    //с аппроксимацией 1-го порядка
    solve_bvp(1, h, n, x, y1);
    
    //с аппроксимацией 2-го порядка
    solve_bvp(2, h, n, x, y2);
    
    //аналитическое решение
    for (int i = 0; i < n; i++) {
        y_analytical[i] = analytical_solution(x[i]);
    }
    
    printf("|  i  |    x    |  Аналит.  |  1-й пор.  |  Погр.1  |  2-й пор.  |  Погр.2  |\n");
    printf("|-----|---------|-----------|------------|----------|------------|----------|\n");
    
    double max_error1 = 0, max_error2 = 0;
    
    for (int i = 0; i < n; i++) {
        double error1 = fabs(y1[i] - y_analytical[i]);
        double error2 = fabs(y2[i] - y_analytical[i]);
        
        if (error1 > max_error1) max_error1 = error1;
        if (error2 > max_error2) max_error2 = error2;
        
        printf("| %3d | %7.3f | %9.4f | %10.4f | %8.4f | %10.4f | %8.4f |\n",
               i, x[i], y_analytical[i], y1[i], error1, y2[i], error2);
    }
    
    printf("\nМаксимальная погрешность:\n");
    printf("Схема 1-го порядка: %.6f\n", max_error1);
    printf("Схема 2-го порядка: %.6f\n", max_error2);
    
    return 0;
}