#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_STAGES 10
#define MAX_METHODS 20
#define MAX_ITER 100
#define TOL 1e-12

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
    double y_exact;
    double error;
} SolutionPoint;

double f(double x, double y) {
    return y*y - y*sin(x) + cos(x);
}

double df_dy(double x, double y) {
    return 2*y - sin(x);
}

double exact_solution(double x) {
    return sin(x);
}

double explicit_rk_step(double x, double y, double h, ButcherTable *table) {
    double k[MAX_STAGES];
    int i, j;
    
    for (i = 0; i < table->s; i++) {
        double x_i = x + table->c[i] * h;
        double y_i = y;
        
        for (j = 0; j < i; j++) {
            y_i += table->a[i][j] * k[j];
        }
        
        k[i] = h * f(x_i, y_i);
    }
    
    double delta_y = 0;
    for (i = 0; i < table->s; i++) {
        delta_y += table->b[i] * k[i];
    }
    
    return y + delta_y;
}

void gauss_solve(int n, double J[MAX_STAGES][MAX_STAGES], double F[MAX_STAGES], double delta_K[MAX_STAGES]) {
    double augmented[MAX_STAGES][MAX_STAGES + 1];
    int i, j, k;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            augmented[i][j] = J[i][j];
        }
        augmented[i][n] = -F[i];
    }
    
    for (k = 0; k < n; k++) {
        double pivot = augmented[k][k];
        for (j = k; j <= n; j++) {
            augmented[k][j] /= pivot;
        }
        
        for (i = 0; i < n; i++) {
            if (i != k) {
                double factor = augmented[i][k];
                for (j = k; j <= n; j++) {
                    augmented[i][j] -= factor * augmented[k][j];
                }
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        delta_K[i] = augmented[i][n];
    }
}

double implicit_rk_step(double x, double y, double h, ButcherTable *table) {
    double K[MAX_STAGES];           
    double K_new[MAX_STAGES];       
    double F[MAX_STAGES];          
    double J[MAX_STAGES][MAX_STAGES]; 
    double delta_K[MAX_STAGES];     
    int i, j, iter;
    
    for (i = 0; i < table->s; i++) {
        K[i] = f(x, y);
    }
    
    for (iter = 0; iter < MAX_ITER; iter++) {
        for (i = 0; i < table->s; i++) {
            double x_i = x + table->c[i] * h;
            double y_i = y;
            
            for (j = 0; j < table->s; j++) {
                y_i += h * table->a[i][j] * K[j];
            }
            
            F[i] = K[i] - f(x_i, y_i);
        }
        
        for (i = 0; i < table->s; i++) {
            double x_i = x + table->c[i] * h;
            double y_i = y;
            
            for (j = 0; j < table->s; j++) {
                y_i += h * table->a[i][j] * K[j];
            }
            
            double dfdy = df_dy(x_i, y_i);
            
            for (j = 0; j < table->s; j++) {
                J[i][j] = (i == j ? 1.0 : 0.0) - h * table->a[i][j] * dfdy;
            }
        }
        
        gauss_solve(table->s, J, F, delta_K);
        
        double max_diff = 0;
        for (i = 0; i < table->s; i++) {
            K_new[i] = K[i] + delta_K[i];
            double diff = fabs(K_new[i] - K[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        for (i = 0; i < table->s; i++) {
            K[i] = K_new[i];
        }
        
        if (max_diff < TOL) break;
    }
    
    double delta_y = 0;
    for (i = 0; i < table->s; i++) {
        delta_y += table->b[i] * K[i];
    }
    
    return y + h * delta_y;
}

double rk_step(double x, double y, double h, ButcherTable *table) {
    if (table->explicit) {
        return explicit_rk_step(x, y, h, table);
    } else {
        return implicit_rk_step(x, y, h, table);
    }
}

SolutionPoint* solve_ode(ButcherTable *table, double x0, double y0, 
                        double x_end, double h, int *n_points) {
    int n = (int)((x_end - x0) / h) + 1;
    *n_points = n;
    
    SolutionPoint *solution = (SolutionPoint*)malloc(n * sizeof(SolutionPoint));
    
    double x = x0;
    double y = y0;
    
    for (int i = 0; i < n; i++) {
        solution[i].x = x;
        solution[i].y_num = y;
        solution[i].y_exact = exact_solution(x);
        solution[i].error = fabs(solution[i].y_num - solution[i].y_exact);
        
        if (i < n - 1) {
            y = rk_step(x, y, h, table);
            x += h;
        }
    }
    
    return solution;
}

void initialize_butcher_tables(ButcherTable *tables) {
    int i, j;
    
    strcpy(tables[0].name, "Явный Эйлер");
    tables[0].s = 1;
    tables[0].explicit = 1;
    tables[0].order = 1;
    tables[0].c[0] = 0.0;
    tables[0].a[0][0] = 0.0;
    tables[0].b[0] = 1.0;
    
    strcpy(tables[1].name, "Неявный Эйлер");
    tables[1].s = 1;
    tables[1].explicit = 0;
    tables[1].order = 1;
    tables[1].c[0] = 1.0;
    tables[1].a[0][0] = 1.0;
    tables[1].b[0] = 1.0;
    
    strcpy(tables[2].name, "Хойн");
    tables[2].s = 2;
    tables[2].explicit = 1;
    tables[2].order = 2;
    tables[2].c[0] = 0.0;      
    tables[2].c[1] = 0.5;     
    tables[2].a[0][0] = 0.0; tables[2].a[0][1] = 0.0;
    tables[2].a[1][0] = 0.5; tables[2].a[1][1] = 0.0; 
    tables[2].b[0] = 0.0;      
    tables[2].b[1] = 1.0;
    
    strcpy(tables[3].name, "Эйлера-Коши");
    tables[3].s = 2;
    tables[3].explicit = 1;
    tables[3].order = 2;
    tables[3].c[0] = 0.0;
    tables[3].c[1] = 1.0;
    tables[3].a[0][0] = 0.0; tables[3].a[0][1] = 0.0;
    tables[3].a[1][0] = 1.0; tables[3].a[1][1] = 0.0;
    tables[3].b[0] = 0.5;
    tables[3].b[1] = 0.5;
    
    strcpy(tables[4].name, "Рунге-Кутта 4");
    tables[4].s = 4;
    tables[4].explicit = 1;
    tables[4].order = 4;
    tables[4].c[0] = 0.0;
    tables[4].c[1] = 0.5;
    tables[4].c[2] = 0.5;
    tables[4].c[3] = 1.0;
    
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            tables[4].a[i][j] = 0.0;
        }
    }
    tables[4].a[1][0] = 0.5;
    tables[4].a[2][1] = 0.5;
    tables[4].a[3][2] = 1.0;
    
    tables[4].b[0] = 1.0/6.0;
    tables[4].b[1] = 1.0/3.0;
    tables[4].b[2] = 1.0/3.0;
    tables[4].b[3] = 1.0/6.0;
    
    strcpy(tables[5].name, "Трапеций");
    tables[5].s = 2;
    tables[5].explicit = 0;
    tables[5].order = 2;
    tables[5].c[0] = 0.0;
    tables[5].c[1] = 1.0;
    tables[5].a[0][0] = 0.0; tables[5].a[0][1] = 0.0;
    tables[5].a[1][0] = 0.5; tables[5].a[1][1] = 0.5;
    tables[5].b[0] = 0.5;
    tables[5].b[1] = 0.5;
}

void save_plot_data(SolutionPoint *solutions[], int num_methods, int n_points, ButcherTable **tables) {
    FILE *file_solutions = fopen("solutions.dat", "w");
    FILE *file_errors = fopen("errors.dat", "w");
    
    if (!file_solutions || !file_errors) {
        printf("Ошибка создания файлов для графиков!\n");
        return;
    }
    
    // Заголовки для файла решений
    fprintf(file_solutions, "# x Точное ");
    for (int i = 0; i < num_methods; i++) {
        fprintf(file_solutions, "%s ", tables[i]->name);
    }
    fprintf(file_solutions, "\n");
    
    // Заголовки для файла ошибок
    fprintf(file_errors, "# x ");
    for (int i = 0; i < num_methods; i++) {
        fprintf(file_errors, "Ошибка_%s ", tables[i]->name);
    }
    fprintf(file_errors, "\n");
    
    // Данные для графиков
    for (int j = 0; j < n_points; j++) {
        // Решения
        fprintf(file_solutions, "%.6f %.6f ", solutions[0][j].x, solutions[0][j].y_exact);
        for (int i = 0; i < num_methods; i++) {
            fprintf(file_solutions, "%.6f ", solutions[i][j].y_num);
        }
        fprintf(file_solutions, "\n");
        
        // Ошибки
        fprintf(file_errors, "%.6f ", solutions[0][j].x);
        for (int i = 0; i < num_methods; i++) {
            fprintf(file_errors, "%.6f ", solutions[i][j].error);
        }
        fprintf(file_errors, "\n");
    }
    
    fclose(file_solutions);
    fclose(file_errors);
    printf("Данные для графиков сохранены в файлы solutions.dat и errors.dat\n");
}

void create_gnuplot_script(int num_methods, ButcherTable **tables) {
    FILE *script = fopen("plot_results.gnu", "w");
    if (!script) {
        printf("Ошибка создания скрипта GNUplot!\n");
        return;
    }
    
    fprintf(script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,10'\n");
    fprintf(script, "set multiplot layout 2,1 title 'Решение ОДУ: y'' = y^2 - y*sin(x) + cos(x)' font 'Arial,12'\n\n");
    
    // Первый график: решения
    fprintf(script, "set title 'Численные решения и точное решение'\n");
    fprintf(script, "set xlabel 'x'\n");
    fprintf(script, "set ylabel 'y(x)'\n");
    fprintf(script, "set grid\n");
    fprintf(script, "set key outside right top\n");
    fprintf(script, "plot 'solutions.dat' using 1:2 with lines linewidth 2 title 'Точное решение (sin(x))'");
    
    // Цвета для разных методов
    char *colors[] = {"red", "blue", "green", "purple", "orange", "brown"};
    char *point_types[] = {"points", "points", "points", "points", "points", "points"};
    
    for (int i = 0; i < num_methods; i++) {
        fprintf(script, ", '' using 1:%d with %s pointtype %d linecolor rgb '%s' title '%s'", 
                i+3, point_types[i], i+1, colors[i], tables[i]->name);
    }
    fprintf(script, "\n\n");
    
    // Второй график: ошибки
    fprintf(script, "set title 'Погрешности численных методов'\n");
    fprintf(script, "set xlabel 'x'\n");
    fprintf(script, "set ylabel 'Погрешность |y_{числ} - y_{точн}|'\n");
    fprintf(script, "set grid\n");
    fprintf(script, "set logscale y\n");
    fprintf(script, "set key outside right top\n");
    fprintf(script, "plot 'errors.dat' using 1:2 with lines linewidth 2 title 'Ошибка %s'", tables[0]->name);
    
    for (int i = 1; i < num_methods; i++) {
        fprintf(script, ", '' using 1:%d with lines linewidth 2 linecolor rgb '%s' title 'Ошибка %s'", 
                i+2, colors[i], tables[i]->name);
    }
    fprintf(script, "\n\n");
    
    fprintf(script, "unset multiplot\n");
    fprintf(script, "set output 'results_comparison.png'\n");
    fprintf(script, "replot\n");
    
    fclose(script);
    printf("Скрипт GNUplot создан: plot_results.gnu\n");
}

void generate_plots() {
    printf("\nГенерация графиков...\n");
    system("gnuplot plot_results.gnu");
    printf("Графики сохранены в файл: results_comparison.png\n");
    printf("Для просмотра графиков выполните: eog results_comparison.png\n");
}

void print_solution(SolutionPoint *solution, int n, ButcherTable *table, double h) {
    printf("\n%s с шагом h = %.3f\n", table->name, h);
    printf("============================================\n");
    printf("   x    Численное   Точное   Погрешность\n");
    printf("--------------------------------------------\n");
    
    double max_error = 0;
    for (int i = 0; i < n; i++) {
        printf("%6.3f   %8.6f  %8.6f  %8.6f\n", 
               solution[i].x, solution[i].y_num, 
               solution[i].y_exact, solution[i].error);
        
        if (solution[i].error > max_error) {
            max_error = solution[i].error;
        }
    }
    printf("--------------------------------------------\n");
    printf("Максимальная погрешность: %.6f\n", max_error);
}

void convergence_analysis(ButcherTable *table, double x0, double y0, double x_end) {
    printf("\nАнализ сходимости %s:\n", table->name);
    printf("h    \tМаксимальная погрешность     Во сколько уменьшилась \n");
    printf("-----------------------------------------------------------\n");
    
    double prev_error = 0;
    double h_values[] = {0.1, 0.05, 0.025, 0.0125};
    
    for (int i = 0; i < 4; i++) {
        double h = h_values[i];
        int n;
        SolutionPoint *sol = solve_ode(table, x0, y0, x_end, h, &n);
        
        double max_error = 0;
        for (int j = 0; j < n; j++) {
            if (sol[j].error > max_error) {
                max_error = sol[j].error;
            }
        }
        
        double ratio = (i > 0) ? prev_error / max_error : 0;
        printf("%.4f        %.10f\t\t   ", h, max_error);
        if (i > 0) printf("%.2f", ratio);
        printf("\n");
        
        prev_error = max_error;
        free(sol);
    }
}

int main() {
    double x0 = 0.0;      
    double y0 = 0.0;      
    double x_end = 4.0;   
    double h = 0.25;     
    
    printf("============================================\n");
    printf("y' = y^2 - y*sin(x) + cos(x)\n");
    printf("y(0) = 0\n");
    printf("y(x) = sin(x)\n");
    printf("Интервал: [0, 4], Шаг: h = %.2f\n", h);
    printf("============================================\n");
    
    ButcherTable tables[6];
    initialize_butcher_tables(tables);
    int num_methods = 6;
    
    // Создаем массив указателей для передачи в функции
    ButcherTable *table_ptrs[6];
    for (int i = 0; i < num_methods; i++) {
        table_ptrs[i] = &tables[i];
    }
    
    // Сохраняем все решения для графиков
    SolutionPoint *solutions[6];
    int n_points;
    
    printf("\nВычисление решений...\n");
    for (int i = 0; i < num_methods; i++) {
        solutions[i] = solve_ode(&tables[i], x0, y0, x_end, h, &n_points);
        print_solution(solutions[i], n_points, &tables[i], h);
    }
    
    // Создаем графики
    save_plot_data(solutions, num_methods, n_points, table_ptrs);
    create_gnuplot_script(num_methods, table_ptrs);
    generate_plots();
    
    // Анализ сходимости
    printf("\n\nАнализ сходимости:\n");
    printf("============================================\n");
    convergence_analysis(&tables[0], x0, y0, x_end);
    convergence_analysis(&tables[2], x0, y0, x_end);
    convergence_analysis(&tables[4], x0, y0, x_end);
    
    // Освобождаем память
    for (int i = 0; i < num_methods; i++) {
        free(solutions[i]);
    }
    
    return 0;
}