#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_STAGES 10
#define MAX_ITER 100
#define TOL 1e-12
#define EPS 0.000001

typedef struct {
    char name[50];
    int s; 
    int explicit; 
    double c[MAX_STAGES];
    double a[MAX_STAGES][MAX_STAGES];
    double b[MAX_STAGES];
    int order;
} ButcherTable;

typedef struct {
    double x;
    double y_num;
    double z_num;
    double y_exact;
    double error;
} SolutionPoint;

//аналитическое решение
double exact_solution(double x) {
    return x * (1 + log(x)) + 0.75 * x * x * x;
}

// y' = z
double f1(double x, double y, double z) {
    return z;
}

// z' = (1/x)z - (1/x²)y + 3x
double f2(double x, double y, double z) {
    return (1.0/x) * z - (1.0/(x*x)) * y + 3.0 * x;
}

//явный метод Рунге-Кутты для системы
void explicit_rk_system_step(double x, double *y, double *z, double h, ButcherTable *table) {
    double k_y[MAX_STAGES], k_z[MAX_STAGES];
    int i, j;
    
    for (i = 0; i < table->s; i++) {
        double x_i = x + table->c[i] * h;
        double y_i = y[0];
        double z_i = z[0];
        
        for (j = 0; j < i; j++) {
            y_i += table->a[i][j] * k_y[j];
            z_i += table->a[i][j] * k_z[j];
        }
        
        k_y[i] = h * f1(x_i, y_i, z_i);
        k_z[i] = h * f2(x_i, y_i, z_i);
    }
    
    double delta_y = 0, delta_z = 0;
    for (i = 0; i < table->s; i++) {
        delta_y += table->b[i] * k_y[i];
        delta_z += table->b[i] * k_z[i];
    }
    
    y[0] += delta_y;
    z[0] += delta_z;
}

void rk_system_step(double x, double *y, double *z, double h, ButcherTable *table) {
    if (table->explicit) {
        explicit_rk_system_step(x, y, z, h, table);
    } else {
        printf("НЕЯВН----------------------------------------------------------");
    }
}

SolutionPoint* solve_system(ButcherTable *table, double x0, double y0, double z0,
                           double x_end, double h, int *n_points) {
    int n = (int)((x_end - x0) / h) + 1;
    *n_points = n;
    
    SolutionPoint *solution = (SolutionPoint*)malloc(n * sizeof(SolutionPoint));
    
    double x = x0;
    double y = y0;
    double z = z0;
    
    for (int i = 0; i < n; i++) {
        solution[i].x = x;
        solution[i].y_num = y;
        solution[i].z_num = z;
        solution[i].y_exact = exact_solution(x);
        solution[i].error = fabs(solution[i].y_num - solution[i].y_exact);
        
        if (i < n - 1) {
            rk_system_step(x, &y, &z, h, table);
            x += h;
        }
    }
    
    return solution;
}

//метод стрельбы - возвращает значение на правом конце
double shooting_function(double alpha, ButcherTable *table) {
    double a = 0.5, b = 4.5;
    double h = 0.20;
    double y0 = 0.247;
    double z0 = alpha;
    
    int n_points;
    SolutionPoint *solution = solve_system(table, a, y0, z0, b, h, &n_points);
    
    double y_final = solution[n_points-1].y_num;
    
    free(solution);
    return y_final;
}

//метод дихотомии
double bisection_method(double target, double alpha_min, double alpha_max, ButcherTable *table) {
    double alpha_mid, f_min, f_mid, f_max;
    int iter = 0;
    
    printf("Поиск параметра alpha методом половинного деления:\n");
    printf("=================================================\n");
    
    f_min = shooting_function(alpha_min, table) - target;
    f_max = shooting_function(alpha_max, table) - target;
    
    printf("F(alpha_min) = F(%.6f) = %.6f\n", alpha_min, f_min);
    printf("F(alpha_max) = F(%.6f) = %.6f\n", alpha_max, f_max);
    
    if (f_min * f_max > 0) {
        printf("ОШИБКА: F(alpha_min) * F(alpha_max) > 0!\n");
        return (alpha_min + alpha_max) / 2;
    }
    
    for (iter = 0; iter < MAX_ITER; iter++) {
        alpha_mid = (alpha_min + alpha_max) / 2;
        f_mid = shooting_function(alpha_mid, table) - target;
        
        printf("Итерация %2d: alpha = %10.6f, y(b) = %10.6f, невязка = %10.6f\n", 
               iter+1, alpha_mid, shooting_function(alpha_mid, table), f_mid);
        
        if (fabs(f_mid) < EPS) {
            printf("Найдено решение с точностью %f\n", EPS);
            break;
        }
        
        if (f_min * f_mid < 0) {
            alpha_max = alpha_mid;
            f_max = f_mid; 
        } else {
            alpha_min = alpha_mid;
            f_min = f_mid;
        }
    }
    
    if (iter == MAX_ITER) {
        printf("Достигнуто максимальное число итераций\n");
    }
    
    return alpha_mid;
}

void initialize_butcher_tables(ButcherTable *tables) {
    strcpy(tables[0].name, "Рунге-Кутта 4 порядка");
    tables[0].s = 4;
    tables[0].explicit = 1;
    tables[0].order = 4;
    
    tables[0].c[0] = 0.0;
    tables[0].c[1] = 0.5;
    tables[0].c[2] = 0.5;
    tables[0].c[3] = 1.0;
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            tables[0].a[i][j] = 0.0;
        }
    }
    
    tables[0].a[1][0] = 0.5;
    tables[0].a[2][1] = 0.5;
    tables[0].a[3][2] = 1.0;
    
    tables[0].b[0] = 1.0/6.0;
    tables[0].b[1] = 1.0/3.0;
    tables[0].b[2] = 1.0/3.0;
    tables[0].b[3] = 1.0/6.0;
}

void print_butcher_table(ButcherTable *table) {
    printf("\n=== %s ===\n", table->name);
    
    printf("c = [");
    for (int i = 0; i < table->s; i++) {
        printf("%.3f", table->c[i]);
        if (i < table->s - 1) printf(", ");
    }
    printf("]\n");
    
    printf("A = \n");
    for (int i = 0; i < table->s; i++) {
        printf("  [");
        for (int j = 0; j < table->s; j++) {
            printf("%6.3f", table->a[i][j]);
            if (j < table->s - 1) printf(", ");
        }
        printf("]\n");
    }
    
    printf("b = [");
    for (int i = 0; i < table->s; i++) {
        printf("%.3f", table->b[i]);
        if (i < table->s - 1) printf(", ");
    }
    printf("]\n");
}

void print_solution(SolutionPoint *solution, int n, ButcherTable *table, double h) {
    printf("\n%s с шагом h = %.3f\n", table->name, h);
    printf("==================================================================\n");
    printf("   x    Численное   Точное   Погрешность     y'\n");
    printf("------------------------------------------------------------------\n");
    
    double max_error = 0;
    for (int i = 0; i < n; i++) {
        printf("%6.2f   %8.6f  %8.6f  %8.6f  %8.6f\n", 
               solution[i].x, solution[i].y_num, 
               solution[i].y_exact, solution[i].error,
               solution[i].z_num);
        
        if (solution[i].error > max_error) {
            max_error = solution[i].error;
        }
    }
    printf("------------------------------------------------------------------\n");
    printf("Максимальная погрешность: %.6f\n", max_error);
}

int main() {
    double a = 0.5, b = 4.5;
    double h = 0.20;
    double y_a = 0.247, y_b_target = 79.612;
    
    printf("РЕШЕНИЕ КРАЕВОЙ ЗАДАЧИ МЕТОДОМ СТРЕЛЬБЫ\n");
    printf("=======================================\n");
    printf("Уравнение: x²y'' - xy' + y = 3x³\n");
    printf("Граничные условия: y(%.1f) = %.3f, y(%.1f) = %.3f\n", a, y_a, b, y_b_target);
    printf("Шаг: h = %.2f\n", h);
    printf("Аналитическое решение: y(x) = x(1 + ln x) + (3/4)x³\n\n");
    
    ButcherTable tables[1];
    initialize_butcher_tables(tables);
    
    printf("Используемый метод:");
    print_butcher_table(&tables[0]);
    
    double alpha_initial = (y_b_target - y_a) / (b - a);
    printf("\nНачальное приближение для alpha: %.6f\n", alpha_initial);
    
    printf("\nТестовый прогон с начальным alpha:\n");
    double test_result = shooting_function(alpha_initial, &tables[0]);
    printf("y(%.1f) = %.6f (цель: %.3f)\n", b, test_result, y_b_target);
    printf("Невязка: %.6f\n\n", test_result - y_b_target);
    
    double alpha_min = alpha_initial - 50.0;
    double alpha_max = alpha_initial + 50.0;
    double alpha_final = bisection_method(y_b_target, alpha_min, alpha_max, &tables[0]);
    
    printf("\nНайденный параметр alpha = %.8f\n", alpha_final);
    
    printf("\nФинальное решение:\n");
    printf("==================\n");
    
    int n_points;
    SolutionPoint *solution = solve_system(&tables[0], a, y_a, alpha_final, b, h, &n_points);
    print_solution(solution, n_points, &tables[0], h);
    
    free(solution);
    return 0;
}