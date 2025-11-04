#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ITER 1000
#define EPS 1e-6

// Исходная система уравнений
void f(double x1, double x2, double *f1, double *f2) {
    *f1 = x1*x1 + x2*x2 - 5*sin(x1) - 9;
    *f2 = 2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3;
}

// Итерирующие функции (метод релаксации с меньшим шагом)
void phi(double x1, double x2, double *phi1, double *phi2) {
    // Используем меньший параметр релаксации для лучшей сходимости
    double lambda1 = 0.02;
    double lambda2 = 0.02;
    
    *phi1 = x1 - lambda1 * (x1*x1 + x2*x2 - 5*sin(x1) - 9);
    *phi2 = x2 - lambda2 * (2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3);
}

// Вычисление частных производных для проверки сходимости
void calculate_jacobian(double x1, double x2, double jacobian[2][2]) {
    double lambda1 = 0.02;
    double lambda2 = 0.02;
    
    // ∂Φ1/∂x1 = 1 - λ1*(2x1 - 5*cos(x1))
    jacobian[0][0] = 1 - lambda1 * (2*x1 - 5*cos(x1));
    
    // ∂Φ1/∂x2 = -λ1*(2x2)
    jacobian[0][1] = -lambda1 * (2*x2);
    
    // ∂Φ2/∂x1 = -λ2*(4x1 + 2x2 - 4 + x2*sin(x1))
    jacobian[1][0] = -lambda2 * (4*x1 + 2*x2 - 4 + x2*sin(x1));
    
    // ∂Φ2/∂x2 = 1 - λ2*(2x1 - 6x2 - cos(x1))
    jacobian[1][1] = 1 - lambda2 * (2*x1 - 6*x2 - cos(x1));
}

// Проверка условия сходимости
int check_convergence(double x1, double x2) {
    double jacobian[2][2];
    calculate_jacobian(x1, x2, jacobian);
    
    // Проверяем условие сходимости: сумма модулей элементов в каждой строке < 1
    double norm1 = fabs(jacobian[0][0]) + fabs(jacobian[0][1]);
    double norm2 = fabs(jacobian[1][0]) + fabs(jacobian[1][1]);
    
    printf("Матрица Якоби:\n");
    printf("∂Φ1/∂x1 = %.6f, ∂Φ1/∂x2 = %.6f\n", jacobian[0][0], jacobian[0][1]);
    printf("∂Φ2/∂x1 = %.6f, ∂Φ2/∂x2 = %.6f\n", jacobian[1][0], jacobian[1][1]);
    printf("Нормы: %.6f, %.6f\n", norm1, norm2);
    
    return (norm1 < 1.0 && norm2 < 1.0);
}

// Метод простой итерации
int simple_iteration(double x1_0, double x2_0, double *x1_result, double *x2_result) {
    double x1_prev = x1_0, x2_prev = x2_0;
    double x1_new, x2_new;
    int iter = 0;
    
    printf("\n=== Начальное приближение: (%.6f, %.6f) ===\n", x1_0, x2_0);
    
    // Проверяем условие сходимости
    if (!check_convergence(x1_0, x2_0)) {
        printf("ВНИМАНИЕ: Условие сходимости не выполняется!\n");
    } else {
        printf("Условие сходимости выполняется.\n");
    }
    
    printf("\n%4s %12s %12s %12s %12s %12s\n", 
           "k", "x1", "x2", "Δx1", "Δx2", "Невязка");
    printf("------------------------------------------------------------\n");
    
    double max_delta;
    
    do {
        // Вычисляем новые значения
        phi(x1_prev, x2_prev, &x1_new, &x2_new);
        
        // Разности
        double dx1 = x1_new - x1_prev;
        double dx2 = x2_new - x2_prev;
        
        // Вычисляем невязку
        double f1, f2;
        f(x1_new, x2_new, &f1, &f2);
        double residual = sqrt(f1*f1 + f2*f2);
        
        printf("%4d %12.6f %12.6f %12.6f %12.6f %12.6f\n", 
               iter, x1_new, x2_new, dx1, dx2, residual);
        
        // Максимальное изменение
        max_delta = fmax(fabs(dx1), fabs(dx2));
        
        // Обновляем значения для следующей итерации
        x1_prev = x1_new;
        x2_prev = x2_new;
        iter++;
        
        // Критерий остановки
        if (iter >= MAX_ITER) {
            printf("Достигнуто максимальное количество итераций!\n");
            break;
        }
        
    } while (max_delta > EPS && iter < 100); // Останавливаемся при малых изменениях
    
    *x1_result = x1_new;
    *x2_result = x2_new;
    
    return iter;
}

int main() {
    printf("РЕШЕНИЕ СИСТЕМЫ НЕЛИНЕЙНЫХ УРАВНЕНИЙ\n");
    printf("Метод простой итерации\n");
    printf("Точность: %.0e\n\n", EPS);
    
    // Более точные начальные приближения
    double initial_approximations[][2] = {
        {-1.3, 1.6},    // Приближение к первому корню
        {-0.9, -2.1},   // Приближение ко второму корню  
        {2.4, 2.5},     // Приближение к третьему корню
        {3.0, -0.9}     // Приближение к четвертому корню
    };
    
    int num_roots = 4;
    double results[4][2];
    int iterations[4];
    
    // Решение для каждого начального приближения
    for (int i = 0; i < num_roots; i++) {
        double x1_0 = initial_approximations[i][0];
        double x2_0 = initial_approximations[i][1];
        
        double x1_result, x2_result;
        int iter_count = simple_iteration(x1_0, x2_0, &x1_result, &x2_result);
        
        results[i][0] = x1_result;
        results[i][1] = x2_result;
        iterations[i] = iter_count;
        
        // Проверка точности найденного корня
        double f1, f2;
        f(x1_result, x2_result, &f1, &f2);
        
        printf("\nКорень %d: x1 = %.6f, x2 = %.6f\n", i+1, x1_result, x2_result);
        printf("Количество итераций: %d\n", iter_count);
        printf("Невязка: f1 = %.2e, f2 = %.2e\n", f1, f2);
        printf("====================================================\n\n");
    }
    
    // Итоговая таблица результатов
    printf("\nИТОГОВЫЕ РЕЗУЛЬТАТЫ:\n");
    printf("-------------------------------------------------------------\n");
    printf("%-8s %-15s %-15s %-12s %-15s\n", "Корень", "x1", "x2", "Итерации", "Невязка");
    printf("-------------------------------------------------------------\n");
    for (int i = 0; i < num_roots; i++) {
        double f1, f2;
        f(results[i][0], results[i][1], &f1, &f2);
        double residual = sqrt(f1*f1 + f2*f2);
        printf("%-8d %-15.6f %-15.6f %-12d %-15.2e\n", 
               i+1, results[i][0], results[i][1], iterations[i], residual);
    }
    
    // Сравнение с ожидаемыми корнями
    printf("\nСРАВНЕНИЕ С ОЖИДАЕМЫМИ КОРНЯМИ:\n");
    printf("-------------------------------------------------------------\n");
    printf("%-8s %-15s %-15s %-15s %-15s\n", "Корень", "Найден x1", "Ожидаемый x1", "Найден x2", "Ожидаемый x2");
    printf("-------------------------------------------------------------\n");
    
    double expected[][2] = {
        {-1.318459, 1.555637},
        {-0.881570, -2.089052},
        {2.423889, 2.532399},
        {2.987711, -0.916494}
    };
    
    for (int i = 0; i < num_roots; i++) {
        printf("%-8d %-15.6f %-15.6f %-15.6f %-15.6f\n", 
               i+1, results[i][0], expected[i][0], results[i][1], expected[i][1]);
    }
    
    return 0;
}