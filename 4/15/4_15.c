#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_STAGES 10
#define MAX_METHODS 20
#define MAX_ITER 100
#define TOL 1e-12
#define EPS 0.0001


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

//частная производная по y
double df_dy(double x, double y) {
    return 2*y - sin(x);
}

//точное решение
double exact_solution(double x) {
    return sin(x);
}

//явный метод Рунге-Кутты
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
    
    // Создаем расширенную матрицу [J | -F]
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            augmented[i][j] = J[i][j];
        }
        augmented[i][n] = -F[i];
    }
    
    // Прямой ход метода Гаусса
    for (k = 0; k < n; k++) {
        // Нормализуем текущую строку
        double pivot = augmented[k][k];
        for (j = k; j <= n; j++) {
            augmented[k][j] /= pivot;
        }
        
        // Исключаем переменную из других строк
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
        delta_K[i] = augmented[i][n];
    }
}

//неявный метод Рунге-Кутты
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
    
    // Метод Ньютона
    for (iter = 0; iter < MAX_ITER; iter++) {
        //вычисляем вектор невязок F
        for (i = 0; i < table->s; i++) {
            double x_i = x + table->c[i] * h;
            double y_i = y;
            
            for (j = 0; j < table->s; j++) {
                y_i += h * table->a[i][j] * K[j];
            }
            
            F[i] = K[i] - f(x_i, y_i);
        }
        
        //вычисляем матрицу Якоби J
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
        
        //решаем систему J·ΔK = -F методом Гаусса
        gauss_solve(table->s, J, F, delta_K);
        
        //обновляем все K 
        double max_diff = 0;
        for (i = 0; i < table->s; i++) {
            K_new[i] = K[i] + delta_K[i];
            double diff = fabs(K_new[i] - K[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        for (i = 0; i < table->s; i++) {
            K[i] = K_new[i];
        }
        
        if (max_diff < TOL) {
            break;
        }
    }
    
    //y_{new} = y + h·∑ b_i·K_i
    double delta_y = 0;
    for (i = 0; i < table->s; i++) {
        delta_y += table->b[i] * K[i];
    }
    
    return y + h * delta_y;
}

//один шаг метода Рунге-Кутты
double rk_step(double x, double y, double h, ButcherTable *table) {
    if (table->explicit) {
        return explicit_rk_step(x, y, h, table);
    } else {
        return implicit_rk_step(x, y, h, table);
    }
}

//решение задачи Коши на всем интервале
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
    
    //явный Эйлер (порядок 1)
    strcpy(tables[0].name, "Явный Эйлер (порядок 1)");
    tables[0].s = 1;
    tables[0].explicit = 1; //явный
    tables[0].order = 1;
    tables[0].c[0] = 0.0;
    tables[0].a[0][0] = 0.0;
    tables[0].b[0] = 1.0;
    
    //неявный Эйлер (порядок 1)
    strcpy(tables[1].name, "Неявный Эйлер (порядок 1)");
    tables[1].s = 1;
    tables[1].explicit = 0;
    tables[1].order = 1;
    tables[1].c[0] = 1.0;
    tables[1].a[0][0] = 1.0;
    tables[1].b[0] = 1.0;
    
    //Хойн (порядок 2)
    strcpy(tables[2].name, "Хойн (порядок 2)");
    tables[2].s = 2;
    tables[2].explicit = 1;
    tables[2].order = 2;
    tables[2].c[0] = 0.0;      
    tables[2].c[1] = 0.5;     
    tables[2].a[0][0] = 0.0; tables[2].a[0][1] = 0.0;
    tables[2].a[1][0] = 0.5; tables[2].a[1][1] = 0.0; 
    tables[2].b[0] = 0.0;      
    tables[2].b[1] = 1.0;      
        
    //Эйлера-Коши (порядок 2)
    strcpy(tables[3].name, "Эйлера-Коши (порядок 2)");
    tables[3].s = 2;
    tables[3].explicit = 1;
    tables[3].order = 2;
    tables[3].c[0] = 0.0;
    tables[3].c[1] = 1.0;
    tables[3].a[0][0] = 0.0; tables[3].a[0][1] = 0.0;
    tables[3].a[1][0] = 1.0; tables[3].a[1][1] = 0.0;
    tables[3].b[0] = 0.5;
    tables[3].b[1] = 0.5;
    
    //классический (порядок 4)
    strcpy(tables[4].name, "Классический (порядок 4)");
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
    
    //метод трапеций неявный (порядок 2)
    strcpy(tables[5].name, "Трапеций неявный (порядок 2)");
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
    
    //инициализация таблиц Бутчера
    ButcherTable tables[6];
    initialize_butcher_tables(tables);
    
    int num_methods = 6;
    
    printf("\nТаблица Бутчера:\n");
    for (int i = 0; i < num_methods; i++) {
        print_butcher_table(&tables[i]);
    }
    
    //решение для каждого метода
    for (int i = 0; i < num_methods; i++) {
        int n_points;
        SolutionPoint *solution = solve_ode(&tables[i], x0, y0, x_end, h, &n_points);
        print_solution(solution, n_points, &tables[i], h);
        free(solution);
    }
    
    printf("\n\nАнализ сходимости:\n");
    printf("============================================\n");
    
    convergence_analysis(&tables[0], x0, y0, x_end); // Явный Эйлер
    convergence_analysis(&tables[2], x0, y0, x_end); // Хойн
    convergence_analysis(&tables[4], x0, y0, x_end); // РК4
    
    return 0;
}