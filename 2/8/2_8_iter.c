#include <stdio.h>
#include <math.h>
#include <float.h>

#define MAX_ITER 1000
#define EPS 0.0001

#define LAMBDA1_X1 0.15
#define LAMBDA1_X2 0.15

#define LAMBDA2_X1 0.03
#define LAMBDA2_X2 0.03

#define LAMBDA3_X1 0.1
#define LAMBDA3_X2 0.1

#define LAMBDA4_X1 0.02
#define LAMBDA4_X2 0.02

double f1(double x1, double x2);
double f2(double x1, double x2);

double phi1_root1(double x1, double x2);
double phi2_root1(double x1, double x2);

double phi1_root2(double x1, double x2);
double phi2_root2(double x1, double x2);

double phi1_root3(double x1, double x2);
double phi2_root3(double x1, double x2);

double phi1_root4(double x1, double x2);
double phi2_root4(double x1, double x2);

int simple_iteration_method(int root_num, double x10, double x20, double *x1_sol, double *x2_sol, int *iterations);

typedef double (*phi_func)(double, double);

int main() {
    printf("\n--- РЕШЕНИЕ СИСТЕМЫ НЕЛИНЕЙНЫХ УРАВНЕНИЙ МЕТОДОМ ПРОСТЫХ ИТЕРАЦИЙ ---\n\n");
    printf("Система уравнений:\n");
    printf("x1^2 + x2^2 - 5*sin(x1) = 9\n");
    printf("2x1^2 + 2x1*x2 - 3x2^2 - 4x1 - x2*cos(x1) + 3 = 0\n\n");

    double initial_approximations[][2] = {
        {-1.3, 1.5},    
        {-0.9, -2.1},   
        {2.4, 2.5},     
        {3.0, -0.9}     
    };

    /*
        x1=-1.318459, x2=1.555637 
        x1=-0.881570, x2=-2.089052 
        x1=2.423889, x2=2.532399 
        x1=2.987711, x2=-0.916494 
    */
    
    int num_approximations = sizeof(initial_approximations) / sizeof(initial_approximations[0]);

    for (int i = 0; i < num_approximations; i++) {
        //берем начальное приближение
        double x10 = initial_approximations[i][0];
        double x20 = initial_approximations[i][1];
        
        printf("Начальное приближение x10=%.1f, x20=%.1f\n", x10, x20);

        double x1_sol, x2_sol;
        int iterations;
        
        int result = simple_iteration_method(i, x10, x20, &x1_sol, &x2_sol, &iterations);

        if (result == 0) {
            printf("  Решение: x1 = %.6f, x2 = %.6f\n", x1_sol, x2_sol);
            printf("  Количество итераций: %d\n\n", iterations);
        } else {
            printf("  Ошибка итерационного метода\n\n");
        }
    }

    return 0;
}

double f1(double x1, double x2) {
    return x1 * x1 + x2 * x2 - 5 * sin(x1) - 9;
}

double f2(double x1, double x2) {
    return 2 * x1 * x1 + 2 * x1 * x2 - 3 * x2 * x2 - 4 * x1 - x2 * cos(x1) + 3;
}

// Корень 1: (-1.318459, 1.555637)
double phi1_root1(double x1, double x2) {
    return x1 + LAMBDA1_X1 * (x1*x1 + x2*x2 - 5*sin(x1) - 9);
}

double phi2_root1(double x1, double x2) {
    return x2 + LAMBDA1_X2 * (2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3);
}

//корень 2: (-0.881570, -2.089052)
double phi1_root2(double x1, double x2) {
    return x1 + LAMBDA2_X1 * (x1*x1 + x2*x2 - 5*sin(x1) - 9);
}

double phi2_root2(double x1, double x2) {
    return x2 + LAMBDA2_X2 * (2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3);
}

//корень 3: (2.423889, 2.532399)
double phi1_root3(double x1, double x2) {
    return x1 + LAMBDA3_X1 * (x1*x1 + x2*x2 - 5*sin(x1) - 9);
}

double phi2_root3(double x1, double x2) {
    return x2 + LAMBDA3_X2 * (2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3);
}

//корень 4: (2.987711, -0.916494)
double phi1_root4(double x1, double x2) {
    return x1 + LAMBDA4_X1 * (x1*x1 + x2*x2 - 5*sin(x1) - 9);
}

double phi2_root4(double x1, double x2) {
    return x2 + LAMBDA4_X2 * (2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3);
}

double phi1_deriv_x1(double x1, double x2, int root_num) {
    double lambda;
    switch(root_num) {
        case 0: lambda = LAMBDA1_X1; break;
        case 1: lambda = LAMBDA2_X1; break;
        case 2: lambda = LAMBDA3_X1; break;
        case 3: lambda = LAMBDA4_X1; break;
        default: lambda = LAMBDA1_X1;
    }
    return 1 + lambda * (2*x1 - 5*cos(x1));
}

double phi1_deriv_x2(double x1, double x2, int root_num) {
    double lambda;
    switch(root_num) {
        case 0: lambda = LAMBDA1_X1; break;
        case 1: lambda = LAMBDA2_X1; break;
        case 2: lambda = LAMBDA3_X1; break;
        case 3: lambda = LAMBDA4_X1; break;
        default: lambda = LAMBDA1_X1;
    }
    return lambda * (2*x2);
}

double phi2_deriv_x1(double x1, double x2, int root_num) {
    double lambda;
    switch(root_num) {
        case 0: lambda = LAMBDA1_X2; break;
        case 1: lambda = LAMBDA2_X2; break;
        case 2: lambda = LAMBDA3_X2; break;
        case 3: lambda = LAMBDA4_X2; break;
        default: lambda = LAMBDA1_X2;
    }
    return lambda * (4*x1 + 2*x2 - 4 + x2*sin(x1));
}

double phi2_deriv_x2(double x1, double x2, int root_num) {
    double lambda;
    switch(root_num) {
        case 0: lambda = LAMBDA1_X2; break;
        case 1: lambda = LAMBDA2_X2; break;
        case 2: lambda = LAMBDA3_X2; break;
        case 3: lambda = LAMBDA4_X2; break;
        default: lambda = LAMBDA1_X2;
    }
    return 1 + lambda * (2*x1 - 6*x2 - cos(x1));
}

//проверка условия сходимости
int check_convergence(double x1, double x2, int root_num, double *max_norm) {
    //вычисляем элементы матрицы Якоби
    double J11 = phi1_deriv_x1(x1, x2, root_num);
    double J12 = phi1_deriv_x2(x1, x2, root_num);
    double J21 = phi2_deriv_x1(x1, x2, root_num);
    double J22 = phi2_deriv_x2(x1, x2, root_num);
    
    if (isnan(J11) || isnan(J12) || isnan(J21) || isnan(J22)) {
        return -1;
    }
    
    double norm1 = fabs(J11) + fabs(J12);
    double norm2 = fabs(J21) + fabs(J22);
    
    //находим максимальную норму
    *max_norm = (norm1 > norm2) ? norm1 : norm2;
    
    // ‖J‖ < 1
    if (*max_norm < 1) {
        return 1;  
    } else {
        return 0;
    }
}

//метод простой итерации 
int simple_iteration_method(int root_num, double x10, double x20, double *x1_sol, double *x2_sol, int *iterations) {
    double x1 = x10;
    double x2 = x20;
    
    // Выбор соответствующих итерирующих функций
    phi_func phi1_func, phi2_func;
    switch(root_num) {
        case 0: 
            phi1_func = phi1_root1;
            phi2_func = phi2_root1;
            break;
        case 1:
            phi1_func = phi1_root2;
            phi2_func = phi2_root2;
            break;
        case 2:
            phi1_func = phi1_root3;
            phi2_func = phi2_root3;
            break;
        case 3:
            phi1_func = phi1_root4;
            phi2_func = phi2_root4;
            break;
        default:
            phi1_func = phi1_root1;
            phi2_func = phi2_root1;
    }
    
    for (int i = 0; i < MAX_ITER; i++) {
        double max_norm;
        int convergence_result = check_convergence(x1, x2, root_num, &max_norm);
        
        if (convergence_result == -1) {
            printf("  Ошибка при проверке сходимости на итерации %d\n", i + 1);
            return -1;
        }
        
        if (convergence_result == 0 && i < 1) {  
            printf("  На итерации %d условие сходимости нарушено (норма = %.3f)\n", i + 1, max_norm);
        }

        double x1_old = x1;
        double x2_old = x2;
        
        //вычисляем новые приближения
        double x1_new = phi1_func(x1, x2);
        double x2_new = phi2_func(x1, x2);
        
        if (isnan(x1_new) || isnan(x2_new)) {
            if (i > 0) {
                x1_new = (x1 + x1_old) / 2;
                x2_new = (x2 + x2_old) / 2;
            } else {
                printf("  Ошибка на итерации %d\n", i + 1);
                return -1;
            }
        }
        
        if (fabs(x1_new - x1) < EPS && fabs(x2_new - x2) < EPS) {
            *x1_sol = x1_new;
            *x2_sol = x2_new;
            *iterations = i + 1;
            return 0;
        }
        
        if (fabs(x1_new) > 100 || fabs(x2_new) > 100) {
            printf("  Расходимость на итерации %d\n", i + 1);
            return -1;
        }
        
        x1 = x1_new;
        x2 = x2_new;
        
        if (root_num == 1 || root_num == 3) {
            if (i % 20 == 0) {
                printf("  Итерация %d: x1=%.6f, x2=%.6f\n", i, x1, x2);
            }
        }
    }
    
    printf("  Метод не сошёлся за %d итераций!\n", MAX_ITER);
    return -1;
}