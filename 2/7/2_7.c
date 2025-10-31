#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS 0.0001
#define MAX_ITER 1000

double f(double x) {
    return 2 * x * x * x - 5 * x * x - sin(3 * x + 2) + 4;
}

double df(double x) {
    return 6 * x * x - 10 * x - 3 * cos(3 * x + 2);
}

double d2f(double x) {
    return 12 * x - 10 + 9 * sin(3 * x + 2);
}

//метод дихотомии
int bisection(double a, double b, double eps, double *root, int *iterations) {
    *iterations = 0;
    if (f(a) * f(b) >= 0) return -1; //условие
    
    while ((b - a) / 2.0 > eps && *iterations < MAX_ITER) {
        double c = (a + b) / 2.0;
        (*iterations)++;
        
        if (f(c) == 0) {
            *root = c;
            return 0;
        }
        else if (f(a) * f(c) < 0) b = c;
        else a = c;
    }
    *root = (a + b) / 2.0;
    return 0;
}

//метод Ньютона
int newton(double a, double b, double eps, double *root, int *iterations) {
    *iterations = 0;
    
    double x;
    double fa = f(a), f2a = d2f(a);
    double fb = f(b), f2b = d2f(b);

    int optimal_a = (fa * f2a > 0) && (fabs(f2a) > 1e-10);
    int optimal_b = (fb * f2b > 0) && (fabs(f2b) > 1e-10);

    if (optimal_a && optimal_b) {
        //обе точки оптимальны - выбираем с большим |f(x)|
        x = (fabs(fa) > fabs(fb)) ? a : b;
    }
    else if (optimal_a) {
        x = a;
    }
    else if (optimal_b) {
        x = b;
    }
    else {
        printf("  ОШИБКА: ни одна точка не удовлетворяет достаточному условию сходимости!\n");
        return -1;
    }
    
    int condition_holds = 1;
    
    double dfx = df(x);
    if (fabs(dfx) < 1e-15) {
        printf("  ВНИМАНИЕ: f'(x₀) слишком мала, сходимость не гарантирована\n");
        condition_holds = 0;
    }
    
    double fx = f(x);
    double d2fx = d2f(x);
    
    //условие
    // |f(x₀)·f''(x₀)| < [f'(x₀)]²
    double left_side = fabs(fx * d2fx);
    double right_side = dfx * dfx;
    
    if (left_side >= right_side) {
        printf("  ВНИМАНИЕ: достаточное условие сходимости не выполняется\n");
        printf("  |f(x₀)·f''(x₀)| = %.6f, [f'(x₀)]² = %.6f\n", left_side, right_side);
        condition_holds = 0;
    }
    
    while (*iterations < MAX_ITER) {
        fx = f(x);
        dfx = df(x);
        
        if (fabs(dfx) < 1e-15) {
            printf("  ОШИБКА: f'(x) слишком близка к 0\n");
            return -1;
        }
        
        double x_new = x - fx / dfx;
        (*iterations)++;
                
        if (fabs(x_new - x) < eps) {
            *root = x_new;
            return 0;
        }
        x = x_new;
    }
    
    return -1;
}

//метод секущих
int secant(double a, double b, double eps, double *root, int *iterations) {
    *iterations = 0;
    
    double x0, x1;
    double fa = f(a), f2a = d2f(a);
    double fb = f(b), f2b = d2f(b);

    int optimal_a = (fa * f2a > 0) && (fabs(f2a) > 1e-10);
    int optimal_b = (fb * f2b > 0) && (fabs(f2b) > 1e-10);

    if (optimal_a && optimal_b) {
        //обе точки оптимальны - выбираем с большим |f(x)|
        x0 = (fabs(fa) > fabs(fb)) ? a : b;
    }
    else if (optimal_a) {
        x0 = a;
    }
    else if (optimal_b) {
        x0 = b;
    }
    else {
        printf("  ОШИБКА: ни одна точка не удовлетворяет достаточному условию сходимости!\n");
        return -1;
    }
    
    // x1 = x0 - f(x0)/f'(x0)
    double dfx0 = df(x0);
    if (fabs(dfx0) < 1e-15) {
        printf("  ОШИБКА: f'(x₀) слишком близка к 0\n");
        return -1;
    }
    x1 = x0 - f(x0) / dfx0;
    
    //условие
    int condition_holds = 1;
    
    for (double test_x = a; test_x <= b; test_x += (b-a)/10.0) {
        double fx = f(test_x);
        double dfx = df(test_x);
        double d2fx = d2f(test_x);
        
        // |f(x)·f''(x)| < [f'(x)]²
        double left_side = fabs(fx * d2fx);
        double right_side = dfx * dfx;
        
        if (left_side >= right_side) {
            printf("  ВНИМАНИЕ: условие сходимости нарушено в точке x=%.3f\n", test_x);
            printf("  |f(x)·f''(x)| = %.6f, [f'(x)]² = %.6f\n", left_side, right_side);
            condition_holds = 0;
            break;
        }
    }
    
    double x_prev = x0;
    double x_curr = x1;
    
    while (*iterations < MAX_ITER) {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);
        
        if (fabs(f_curr - f_prev) < 1e-15) {
            printf("  ОШИБКА\n");
            return -1;
        }
        
        double x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
        (*iterations)++;
                
        if (fabs(x_next - x_curr) < eps) {
            *root = x_next;
            return 0;
        }
        
        x_prev = x_curr;
        x_curr = x_next;
    }
    
    return -1;
}

//метод хорд
int chord(double a, double b, double eps, double *root, int *iterations) {
    *iterations = 0;
    if (f(a) * f(b) >= 0) {
        printf("  ОШИБКА: на интервале [%.2f, %.2f] нет корня\n", a, b);
        return -1;
    }

    double z, x;
    double fa = f(a), f2a = d2f(a);
    double fb = f(b), f2b = d2f(b);

    int optimal_a = (fa * f2a > 0) && (fabs(f2a) > 1e-10);
    int optimal_b = (fb * f2b > 0) && (fabs(f2b) > 1e-10);

    if (optimal_a && optimal_b) {
        // Обе точки оптимальны
        z = (fabs(fa) > fabs(fb)) ? a : b;  // z - с большим |f(x)|
        x = (z == a) ? b : a;               // x - другая точка
    }
    else if (optimal_a) {
        z = a;  // a оптимальна для z
        x = b;  // x - другая точка
    }
    else if (optimal_b) {
        z = b;  // b оптимальна для z  
        x = a;  // x - другая точка
    }
    else {
        printf("  ОШИБКА: ни одна точка не удовлетворяет достаточному условию сходимости!\n");
        return -1;
    }
    
    double fz = f(z);  

    int condition_holds = 1;
    
    //условие
    for (double test_x = a; test_x <= b; test_x += (b-a)/10.0) {
        double fx = f(test_x);
        double dfx = df(test_x);
        double d2fx = d2f(test_x);
        
        // |f(x)·f''(x)| < [f'(x)]²
        double left_side = fabs(fx * d2fx);
        double right_side = dfx * dfx;
        
        if (left_side >= right_side) {
            printf("  ВНИМАНИЕ: условие сходимости нарушено в точке x=%.3f\n", test_x);
            condition_holds = 0;
            break;
        }
    }
    
    if (!condition_holds) {
        printf("  ВНИМАНИЕ: достаточное условие сходимости не выполняется на интервале\n");
    }
    
    while (*iterations < MAX_ITER) {
        double fx = f(x);  
        
        double x_new = x - fx * (z - x) / (fz - fx);
        (*iterations)++;
                
        if (fabs(x_new - x) < eps) {
            *root = x_new;
            return 0;
        }
        
        x = x_new; 
    }
    
    return -1;
}

//метод простой итерации
int simple_iteration(double a, double b, double eps, double *root, int *iterations) {
    *iterations = 0;
    
    double x = (a + b) / 2.0;
    
    double lambda;
    double df_x0 = df(x);
    
    if (fabs(df_x0) > 1e-10) {
        lambda = -0.5 / df_x0;  
    } else {
        lambda = -0.05;  
    }
    
    if (fabs(lambda) >= 0.5) {
        lambda = (lambda > 0) ? 0.49 : -0.49;
    }
    
    //условие
    int condition_holds = 1;
    for (double test_x = a; test_x <= b; test_x += (b-a)/10.0) {
        double phi_prime = 1 + lambda * df(test_x);
        if (fabs(phi_prime) >= 1.0) {
            condition_holds = 0;
            break;
        }
    }
    
    if (!condition_holds) {
        printf("  ОШИБКА: условие |φ'(x)| < 1 не выполняется на интервале!\n");
        return -1;
    }
    
    double x_prev = x;
    double x_curr = x_prev + lambda * f(x_prev);
    (*iterations)++;
    
    while (*iterations < MAX_ITER) {
        if (fabs(x_curr - x_prev) < eps) {
            *root = x_curr;
            return 0;
        }
        
        x_prev = x_curr;
        x_curr = x_prev + lambda * f(x_prev);
        (*iterations)++;
    }
    
    return -1;
}

//поиск всех интервалов, содержащих корни
void find_root_intervals(double x_min, double x_max, double step, 
                        double intervals[][2], int *interval_count) {
    *interval_count = 0;
    double x_left = x_min;
    double f_left = f(x_left);
    
    for (double x = x_min + step; x <= x_max && *interval_count < 10; x += step) {
        double f_current = f(x);
        
        if (f_left * f_current <= 0) {
            // Найден интервал с корнем
            intervals[*interval_count][0] = x_left;
            intervals[*interval_count][1] = x;
            (*interval_count)++;
        }
        x_left = x;
        f_left = f_current;
    }
}

//генерация данных функции для построения графика
void generate_function_data(const char* filename, double x_min, double x_max, int points) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Не удалось создать файл данных функции");
        return;
    }
    double dx = (x_max - x_min) / (points - 1);
    for (int i = 0; i < points; i++) {
        double x = x_min + i * dx;
        fprintf(fp, "%.6f %.6f\n", x, f(x));
    }
    fclose(fp);
}

//сохранение корней в файл
void save_roots(const char* filename, double roots[], int count) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Не удалось создать файл корней");
        return;
    }
    for (int i = 0; i < count; i++) {
        fprintf(fp, "%.6f %.6f\n", roots[i], f(roots[i]));
    }
    fclose(fp);
}

//построение графика через gnuplot
void plot_with_gnuplot(double roots[], int root_count) {
    printf("\nПостроение графика...\n");
    
    //генерируем данные функции
    generate_function_data("function_data.txt", -2.0, 3.0, 500);
    save_roots("roots.txt", roots, root_count);

    //создаём скрипт gnuplot
    FILE *gp = fopen("plot.gp", "w");
    if (!gp) {
        perror("Не удалось создать скрипт gnuplot");
        return;
    }

    fprintf(gp, "set terminal pngcairo size 1000,600 enhanced font 'Arial,12'\n");
    fprintf(gp, "set output 'equation_plot.png'\n");
    fprintf(gp, "set title 'График функции f(x) = 2x^{3} - 5x^{2} - sin(3x+2) + 4'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'f(x)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top left\n");
    fprintf(gp, "set zeroaxis lt -1 lc 'black'\n");
    fprintf(gp, "plot 'function_data.txt' with lines lw 2 lc rgb 'blue' title 'f(x)', \\\n");
    fprintf(gp, "     'roots.txt' with points pt 7 ps 1.5 lc rgb 'red' title 'Найденные корни (%d шт.)'\n", root_count);
    
    fclose(gp);

    //запускаем gnuplot
    int result = system("gnuplot plot.gp");
    if (result == 0) {
        printf("График сохранен в файл 'equation_plot.png'\n");
        
        // Показываем график
        printf("Открываю график...\n");
        system("display equation_plot.png 2>/dev/null || xdg-open equation_plot.png 2>/dev/null || open equation_plot.png 2>/dev/null");
    } else {
        printf("Ошибка: не удалось запустить gnuplot. Убедитесь, что он установлен.\n");
    }
}

//решение уравнения одним методом для данного интервала
void solve_for_interval(double a, double b, double x0, const char* interval_name) {
    printf("\n=== %s [%.2f, %.2f] ===\n", interval_name, a, b);
    
    double root;
    int iterations;
    
    //проверка наличия корня в интервале
    if (f(a) * f(b) > 0) {
        printf("Нет корня в данном интервале (f(%.2f)=%.2f, f(%.2f)=%.2f)\n", 
               a, f(a), b, f(b));
        return;
    }
    
    // Метод дихотомии
    if (bisection(a, b, EPS, &root, &iterations) == 0) {
        printf("Метод дихотомии:    x = %.6f, f(x) = %.2e, итераций: %d\n", root, f(root), iterations);
    }
    
    // Метод Ньютона
    if (newton(a, b, EPS, &root, &iterations) == 0) {
        printf("Метод Ньютона:      x = %.6f, f(x) = %.2e, итераций: %d\n", root, f(root), iterations);
    }
    
    // Метод секущих
    if (secant(a, b, EPS, &root, &iterations) == 0) {
        printf("Метод секущих:      x = %.6f, f(x) = %.2e, итераций: %d\n", root, f(root), iterations);
    }
    
    // Метод хорд
    if (chord(a, b, EPS, &root, &iterations) == 0) {
        printf("Метод хорд:         x = %.6f, f(x) = %.2e, итераций: %d\n", root, f(root), iterations);
    }
    
    // Метод простой итерации
    if (simple_iteration(a, b, EPS, &root, &iterations) == 0) {
        printf("Метод простой итер.: x = %.6f, f(x) = %.2e, итераций: %d\n", root, f(root), iterations);
    }
}

int main() {
    printf("===============================================\n");
    printf("Решение уравнения: 2x^3 - 5x^2 - sin(3x+2) + 4 = 0\n");
    printf("Точность: %.0e\n", EPS);
    printf("===============================================\n\n");

    //поиск интервалов с корнями
    double intervals[10][2];
    int interval_count;
    
    printf("Поиск интервалов с корнями...\n");
    find_root_intervals(-2.0, 3.0, 0.1, intervals, &interval_count);
    
    printf("Найдено интервалов с корнями: %d\n", interval_count);
    for (int i = 0; i < interval_count; i++) {
        printf("Интервал %d: [%.2f, %.2f], f(a)=%.2f, f(b)=%.2f\n", 
               i+1, intervals[i][0], intervals[i][1], 
               f(intervals[i][0]), f(intervals[i][1]));
    }
    printf("\n");

    double roots[10];
    int root_count = 0;
    
    //решение для каждого найденного интервала
    for (int i = 0; i < interval_count; i++) {
        double a = intervals[i][0];
        double b = intervals[i][1];
        double x0 = (a + b) / 2.0; //середина интервала
        
        char interval_name[50];
        sprintf(interval_name, "КОРЕНЬ %d", i+1);
        
        //решение всеми методами
        solve_for_interval(a, b, x0, interval_name);
        
        //сохраняем найденный корень (берем результат метода дихотомии как наиболее надежный)
        double root;
        int iter;
        if (bisection(a, b, EPS, &root, &iter) == 0) {
            roots[root_count++] = root;
        }
    }

    return 0;
}