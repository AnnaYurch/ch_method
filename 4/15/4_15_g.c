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
} SolutionPoint;

// Forward declarations
double f(double x, double y);
double df_dy(double x, double y);
double exact_solution(double x);
double explicit_rk_step(double x, double y, double h, ButcherTable *table);
void gauss_solve(int n, double J[MAX_STAGES][MAX_STAGES], double F[MAX_STAGES], double delta_K[MAX_STAGES]);
double implicit_rk_step(double x, double y, double h, ButcherTable *table);
double rk_step(double x, double y, double h, ButcherTable *table);
SolutionPoint* solve_ode(ButcherTable *table, double x0, double y0, 
                        double x_end, double h, int *n_points);
void initialize_butcher_tables(ButcherTable *tables);
void save_plot_data(SolutionPoint *solutions[], int num_methods, int n_points, ButcherTable **tables);
void create_gnuplot_script(int num_methods, ButcherTable **tables);
void generate_plots();
void print_solution(SolutionPoint *solution, int n, ButcherTable *table, double h);

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
    if (solution == NULL) {
        printf("Ошибка выделения памяти!\n");
        exit(1);
    }
    
    double x = x0;
    double y = y0;
    
    for (int i = 0; i < n; i++) {
        solution[i].x = x;
        solution[i].y_num = y;
        solution[i].y_exact = exact_solution(x);
        
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
    
    if (!file_solutions) {
        printf("Ошибка создания файла solutions.dat!\n");
        return;
    }
    
    // Заголовки для файла решений
    fprintf(file_solutions, "# x exact");
    for (int i = 0; i < num_methods; i++) {
        fprintf(file_solutions, " %s", tables[i]->name);
    }
    fprintf(file_solutions, "\n");
    
    // Данные для графиков
    for (int j = 0; j < n_points; j++) {
        fprintf(file_solutions, "%.6f %.6f", solutions[0][j].x, solutions[0][j].y_exact);
        for (int i = 0; i < num_methods; i++) {
            fprintf(file_solutions, " %.6f", solutions[i][j].y_num);
        }
        fprintf(file_solutions, "\n");
    }
    
    fclose(file_solutions);
    printf("Данные сохранены в solutions.dat: %d точек от x=%.2f до x=%.2f\n", 
           n_points, solutions[0][0].x, solutions[0][n_points-1].x);
}

void create_gnuplot_script(int num_methods, ButcherTable **tables) {
    FILE *script = fopen("plot_results.gnu", "w");
    if (!script) {
        printf("Ошибка создания скрипта GNUplot!\n");
        return;
    }
    
    fprintf(script, "set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'\n");
    fprintf(script, "set output 'ode_solution.png'\n");
    fprintf(script, "set title 'Решение ОДУ: y'' = y^2 - y*sin(x) + cos(x)' font 'Verdana,14'\n");
    fprintf(script, "set xlabel 'x' font 'Verdana,12'\n");
    fprintf(script, "set ylabel 'y(x)' font 'Verdana,12'\n");
    fprintf(script, "set xrange [0:4]\n");
    fprintf(script, "set yrange [-1.2:1.2]\n");
    fprintf(script, "set grid\n");
    fprintf(script, "set key outside right top vertical box\n");
    fprintf(script, "set key spacing 1.5\n\n");
    
    // Точное решение
    fprintf(script, "plot 'solutions.dat' using 1:2 with lines linewidth 3 linecolor rgb 'black' title 'Точное решение (sin(x))'");
    
    // Численные методы
    char *colors[] = {"red", "blue", "green", "purple", "orange", "brown"};
    char *styles[] = {"points pointtype 7 pointsize 1.2", 
                     "points pointtype 9 pointsize 1.2",
                     "points pointtype 11 pointsize 1.2",
                     "points pointtype 13 pointsize 1.2", 
                     "points pointtype 15 pointsize 1.2",
                     "points pointtype 17 pointsize 1.2"};
    
    for (int i = 0; i < num_methods; i++) {
        fprintf(script, ", \\\n     '' using 1:%d with %s linecolor rgb '%s' title '%s'", 
                i+3, styles[i], colors[i], tables[i]->name);
    }
    fprintf(script, "\n");
    
    fclose(script);
    printf("Скрипт GNUplot создан: plot_results.gnu\n");
}

void generate_plots() {
    printf("\nЗапуск GNUplot для создания графика...\n");
    int result = system("gnuplot plot_results.gnu");
    if (result == 0) {
        printf("График успешно создан: ode_solution.png\n");
        
        // Попробуем открыть график разными способами
        printf("Попытка открыть график...\n");
        int opened = system("xdg-open ode_solution.png 2>/dev/null");
        if (opened != 0) {
            opened = system("eog ode_solution.png 2>/dev/null");
        }
        if (opened != 0) {
            opened = system("display ode_solution.png 2>/dev/null");
        }
        if (opened != 0) {
            printf("График сохранен как 'ode_solution.png'. Откройте его вручную.\n");
        }
    } else {
        printf("Ошибка при создании графика!\n");
        printf("Проверьте установлен ли GNUplot: sudo apt install gnuplot\n");
    }
}

void print_solution(SolutionPoint *solution, int n, ButcherTable *table, double h) {
    printf("\n%s (шаг h = %.3f)\n", table->name, h);
    printf(" x     Численное   Точное\n");
    printf("--------------------------\n");
    
    // Показываем только каждую 4-ю точку для компактности
    for (int i = 0; i < n; i += 4) {
        printf("%.3f  %8.6f  %8.6f\n", 
               solution[i].x, solution[i].y_num, solution[i].y_exact);
    }
    if (n > 1 && (n-1) % 4 != 0) {
        printf("%.3f  %8.6f  %8.6f\n", 
               solution[n-1].x, solution[n-1].y_num, solution[n-1].y_exact);
    }
    printf("--------------------------\n");
}

int main() {
    double x0 = 0.0;      
    double y0 = 0.0;      
    double x_end = 4.0;   
    double h = 0.1;
    
    printf("============================================\n");
    printf("Решение ОДУ: y' = y^2 - y*sin(x) + cos(x)\n");
    printf("Начальные условия: y(0) = 0\n");
    printf("Точное решение: y(x) = sin(x)\n");
    printf("Интервал: [%.1f, %.1f], Шаг: h = %.2f\n", x0, x_end, h);
    printf("============================================\n");
    
    ButcherTable tables[6];
    initialize_butcher_tables(tables);
    int num_methods = 6;
    
    // Создаем массив указателей
    ButcherTable *table_ptrs[6];
    for (int i = 0; i < num_methods; i++) {
        table_ptrs[i] = &tables[i];
    }
    
    // Вычисляем решения
    SolutionPoint *solutions[6];
    int n_points;
    
    printf("\nВычисление решений методами Рунге-Кутты...\n");
    for (int i = 0; i < num_methods; i++) {
        solutions[i] = solve_ode(&tables[i], x0, y0, x_end, h, &n_points);
        print_solution(solutions[i], n_points, &tables[i], h);
    }
    
    // Проверяем данные
    printf("\nПроверка данных:\n");
    for (int i = 0; i < num_methods; i++) {
        printf("%s: x ∈ [%.3f, %.3f], точек: %d\n", 
               tables[i].name, solutions[i][0].x, solutions[i][n_points-1].x, n_points);
    }
    
    // Создаем графики
    printf("\nСоздание графиков...\n");
    save_plot_data(solutions, num_methods, n_points, table_ptrs);
    create_gnuplot_script(num_methods, table_ptrs);
    generate_plots();
    
    // Проверяем существование файлов
    printf("\nПроверка созданных файлов:\n");
    system("ls -la solutions.dat plot_results.gnu ode_solution.png 2>/dev/null");
    
    // Освобождаем память
    for (int i = 0; i < num_methods; i++) {
        free(solutions[i]);
    }
    
    printf("\nПрограмма завершена.\n");
    return 0;
}