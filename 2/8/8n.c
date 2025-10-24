#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS 0.0001
#define MAX_ITER 1000

int UL_decomposition(int n, double A[n][n], double LU[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            LU[i][j] = A[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(LU[i][i]) < 1e-12) {
                printf("Ошибка: нулевой диагональный элемент U[%d][%d]\n", i, i);
                return 0;
            }
            LU[j][i] /= LU[i][i]; 
            
            for (int k = i + 1; k < n; k++) {
                LU[j][k] -= LU[j][i] * LU[i][k];
            }
        }
    }
    return 1;
}

// Ly = b
void solve_L(int n, double LU[n][n], double b[n], double y[n]) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];
        }
    }
}

// Ux = y
void solve_U(int n, double LU[n][n], double y[n], double x[n]) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= LU[i][i];
    }
}

int inverse_matrix(int n, double A[n][n], double invA[n][n]) {
    double LU[n][n];
    double y[n], x[n];
    double e[n];
    
    if (!UL_decomposition(n, A, LU)) {
        return 0;
    }
    
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            e[i] = (i == j) ? 1.0 : 0.0;
        }
        
        solve_L(n, LU, e, y);
        solve_U(n, LU, y, x);
        
        for (int i = 0; i < n; i++) {
            invA[i][j] = x[i];
        }
    }
    return 1;
}

void f(double x1, double x2, double *f1, double *f2) {
    *f1 = x1*x1 + x2*x2 - 5*sin(x1) - 9;
    *f2 = 2*x1*x1 + 2*x1*x2 - 3*x2*x2 - 4*x1 - x2*cos(x1) + 3;
}

void jacobian(double x1, double x2, double J[2][2]) {
    J[0][0] = 2*x1 - 5*cos(x1);
    J[0][1] = 2*x2;
    J[1][0] = 4*x1 + 2*x2 - 4 + x2*sin(x1);
    J[1][1] = 2*x1 - 6*x2 - cos(x1);
}

int newton_method(double *x1, double *x2) {
    double x1_old, x2_old;
    double f1, f2;
    double J[2][2], J_inv[2][2];
    double dx[2];
    int iter = 0;
    double dx_norm = 1.0;
    
    printf("Метод Ньютона (с LU-разложением):\n");
    printf("k\tx1\t\tx2\t\tf1\t\tf2\n");
    
    do {
        x1_old = *x1;
        x2_old = *x2;
        
        f(*x1, *x2, &f1, &f2);
        
        if (fabs(f1) < EPS && fabs(f2) < EPS) break;
        
        jacobian(*x1, *x2, J);
        
        if (!inverse_matrix(2, J, J_inv)) {
            dx[0] = -0.01 * f1;
            dx[1] = -0.01 * f2;
        } else {
            dx[0] = -(J_inv[0][0] * f1 + J_inv[0][1] * f2);
            dx[1] = -(J_inv[1][0] * f1 + J_inv[1][1] * f2);
        }
        
        *x1 += dx[0];
        *x2 += dx[1];
        
        dx_norm = sqrt((*x1 - x1_old)*(*x1 - x1_old) + 
                      (*x2 - x2_old)*(*x2 - x2_old));
        
        printf("%d\t%.6f\t%.6f\t%.2e\t%.2e\n", iter, *x1, *x2, f1, f2);
        
        if (++iter >= MAX_ITER) break;
        
    } while (dx_norm > EPS);
    
    return iter;
}

// НОВАЯ функция проверки сходимости с подбором λ
int check_simple_iteration_convergence(double x1, double x2, double lambda1, double lambda2) {
    double J[2][2];
    jacobian(x1, x2, J);
    
    // Матрица производных итерирующих функций: Φi = xi + λi·fi
    double Phi_J[2][2];
    Phi_J[0][0] = 1 + lambda1 * J[0][0];  // dΦ1/dx1
    Phi_J[0][1] = lambda1 * J[0][1];      // dΦ1/dx2
    Phi_J[1][0] = lambda2 * J[1][0];      // dΦ2/dx1  
    Phi_J[1][1] = 1 + lambda2 * J[1][1];  // dΦ2/dx2
    
    // Норма матрицы (кубическая)
    double norm1 = fabs(Phi_J[0][0]) + fabs(Phi_J[0][1]);
    double norm2 = fabs(Phi_J[1][0]) + fabs(Phi_J[1][1]);
    double norm = (norm1 > norm2) ? norm1 : norm2;
        
    if (norm < 1.0) {
        printf("Условие сходимости ВЫПОЛНЕНО (||J_Φ|| = %.6f < 1) при λ=(%.4f,%.4f)\n", norm, lambda1, lambda2);
        return 1;
    } else {
        printf("Условие сходимости НЕ ВЫПОЛНЕНО (||J_Φ|| = %.6f >= 1) при λ=(%.4f,%.4f)\n", norm, lambda1, lambda2);
        return 0;
    }
}

// НОВАЯ функция автоматического подбора λ
void find_optimal_lambda(double x1, double x2, double *lambda1, double *lambda2) {
    double best_lambda1 = 0.001, best_lambda2 = 0.001;
    double best_norm = 1000.0;
    
    // Пробуем разные значения λ
    double lambda_values[] = {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, -0.001, -0.01};
    int num_lambdas = sizeof(lambda_values) / sizeof(lambda_values[0]);
    
    for (int i = 0; i < num_lambdas; i++) {
        for (int j = 0; j < num_lambdas; j++) {
            double lambda1_test = lambda_values[i];
            double lambda2_test = lambda_values[j];
            
            double J[2][2];
            jacobian(x1, x2, J);
            
            double Phi_J[2][2];
            Phi_J[0][0] = 1 + lambda1_test * J[0][0];
            Phi_J[0][1] = lambda1_test * J[0][1];
            Phi_J[1][0] = lambda2_test * J[1][0];
            Phi_J[1][1] = 1 + lambda2_test * J[1][1];
            
            double norm1 = fabs(Phi_J[0][0]) + fabs(Phi_J[0][1]);
            double norm2 = fabs(Phi_J[1][0]) + fabs(Phi_J[1][1]);
            double norm = (norm1 > norm2) ? norm1 : norm2;
            
            if (norm < best_norm && norm < 0.9) {
                best_norm = norm;
                best_lambda1 = lambda1_test;
                best_lambda2 = lambda2_test;
            }
        }
    }
    
    *lambda1 = best_lambda1;
    *lambda2 = best_lambda2;
    printf("Оптимальные параметры: λ1=%.6f, λ2=%.6f (||J_Φ||=%.6f)\n", best_lambda1, best_lambda2, best_norm);
}

// НОВАЯ итерационная функция с разными методами
void phi_simple_theory(double x1, double x2, double *x1_new, double *x2_new, int method_type) {
    double f1, f2;
    f(x1, x2, &f1, &f2);
    
    double lambda1, lambda2;
    
    if (method_type == 1) {
        // Для корней в правой полуплоскости (положительные x1)
        lambda1 = -0.01;
        lambda2 = -0.01;
    }
    else if (method_type == 2) {
        // Для корней в левой полуплоскости (отрицательные x1)
        lambda1 = 0.01;  
        lambda2 = 0.01;
    }
    else if (method_type == 3) {
        // Смешанный вариант
        lambda1 = -0.005;
        lambda2 = 0.005;
    }
    else {
        // Автоматический подбор
        find_optimal_lambda(x1, x2, &lambda1, &lambda2);
    }
    
    *x1_new = x1 + lambda1 * f1;
    *x2_new = x2 + lambda2 * f2;
}

// НОВАЯ улучшенная функция метода простых итераций
int simple_iteration_adaptive(double *x1, double *x2) {
    double x1_old, x2_old;
    int iter = 0;
    double dx_norm = 1.0;
    
    printf("Метод простой итерации:\n");
    printf("Начальное приближение: x1 = %f, x2 = %f\n", *x1, *x2);
    
    // Определяем тип метода по начальному приближению
    int method_type;
    if (*x1 > 0) {
        method_type = 1;  // Для положительных x1
        printf("Используется метод для положительных x1\n");
    } else {
        method_type = 2;  // Для отрицательных x1  
        printf("Используется метод для отрицательных x1\n");
    }
    
    // Проверяем сходимость
    double lambda1, lambda2;
    if (method_type == 1) {
        lambda1 = -0.01; lambda2 = -0.01;
    } else {
        lambda1 = 0.01; lambda2 = 0.01;
    }
    check_simple_iteration_convergence(*x1, *x2, lambda1, lambda2);
    
    for (iter = 0; iter < MAX_ITER; iter++) {
        x1_old = *x1;
        x2_old = *x2;
        
        // Вычисляем x^(k+1) = Φ(x^(k))
        double x1_new, x2_new;
        phi_simple_theory(*x1, *x2, &x1_new, &x2_new, method_type);
        
        if (!isfinite(x1_new) || !isfinite(x2_new)) {
            printf("Расходимость на итерации %d!\n", iter);
            return -1;
        }
        
        *x1 = x1_new;
        *x2 = x2_new;
        
        double f1, f2;
        f(*x1, *x2, &f1, &f2);
        
        dx_norm = sqrt((*x1 - x1_old)*(*x1 - x1_old) + 
                      (*x2 - x2_old)*(*x2 - x2_old));
        
        double residual = sqrt(f1*f1 + f2*f2);
        
        if (iter % 100 == 0 && iter > 0) {
            printf("Итерация %d: x1 = %.8f, x2 = %.8f, ошибка = %e, невязка = %e\n", 
                   iter, *x1, *x2, dx_norm, residual);
        }
        
        if (dx_norm < EPS && residual < EPS) {
            printf("Сходимость достигнута на итерации %d\n", iter);
            printf("Корень: x1 = %.6f, x2 = %.6f\n", *x1, *x2);
            printf("Невязка: f1 = %e, f2 = %e\n", f1, f2);
            return iter;
        }
        
        // Адаптация метода если медленно сходится
        if (iter > 50 && dx_norm > 0.1) {
            printf("Медленная сходимость, пробуем другой метод...\n");
            method_type = 3;  // Переключаемся на смешанный метод
        }
    }
    
    printf("Достигнут предел итераций (%d)\n", MAX_ITER);
    return iter;
}

void phi_seidel_theory(double x1, double x2, double *x1_new, double *x2_new) {
    double f1, f2;
    double lambda1 = 0.001, lambda2 = 0.001;
    
    f(x1, x2, &f1, &f2);
    *x1_new = x1 + lambda1 * f1;
    
    f(*x1_new, x2, &f1, &f2); //подставляем вычесленное x1_new
    *x2_new = x2 + lambda2 * f2;
}

int seidel_method_adaptive(double *x1, double *x2) {
    double x1_old, x2_old;
    double x1_new, x2_new;
    int iter = 0;
    double dx_norm = 1.0;
    
    printf("Метод Зейделя:\n");
    //printf("k\tx1\t\tx2\t\tf1\t\tf2\t\t||dx||\n");
    
    do {
        x1_old = *x1;
        x2_old = *x2;
        
        phi_seidel_theory(*x1, *x2, &x1_new, &x2_new);
        
        if (!isfinite(x1_new) || !isfinite(x2_new)) {
            printf("Метод расходится!\n");
            return -1;
        }
        
        *x1 = x1_new;
        *x2 = x2_new;
        
        double f1, f2;
        f(*x1, *x2, &f1, &f2);
        
        dx_norm = sqrt((*x1 - x1_old)*(*x1 - x1_old) + 
                      (*x2 - x2_old)*(*x2 - x2_old));
        
        //printf("%d\t%.6f\t%.6f\t%.2e\t%.2e\t%.2e\n", 
            //   iter, *x1, *x2, f1, f2, dx_norm);
        
        if (++iter >= MAX_ITER) {
            printf("Достигнут предел итераций\n");
            break;
        }
        
    } while (dx_norm > EPS);
    
    return iter;
}

//функция для поиска точек пересечения графиков эквивалентных функций
//функция для поиска точек пересечения графиков эквивалентных функций
void find_fixed_points(double lambda1, double lambda2) {
    printf("Поиск точек пересечения графиков эквивалентных функций:\n");

    double test_points[][2] = {
        {2.5, 2.5}, {-1.0, -2.0}, {-1.5, 1.5},
        {0.0, 0.0}, {1.0, 1.0}, {-1.0, 1.0},
        {3.0, 3.0}, {-2.0, -3.0}, {-2.0, 2.0}
    };
    
    for (int i = 0; i < 9; i++) {
        double x1 = test_points[i][0];
        double x2 = test_points[i][1];
        double phi1, phi2;
        
        // Определяем method_type по x1 (как в основной функции)
        int method_type;
        if (x1 > 0) {
            method_type = 1;  // Для положительных x1
        } else {
            method_type = 2;  // Для отрицательных x1
        }
        
        phi_simple_theory(x1, x2, &phi1, &phi2, method_type);
        
        double f1, f2;
        f(x1, x2, &f1, &f2);
        
        printf("Точка (%5.1f, %5.1f): ", x1, x2);
        printf("f1=%.3f, f2=%.3f", f1, f2);
        
        //проверка на близость к неподвижной точке
        if (fabs(phi1 - x1) < 0.1 && fabs(phi2 - x2) < 0.1) {
            printf(" - БЛИЗКО К ПЕРЕСЕЧЕНИЮ\n");
        } else {
            printf(" - далеко\n");
        }
    }
}

//==============

int main() {

    find_fixed_points(0.001, 0.001);

    /*
    double solutions[][2] = {
        {2.0, 2.5},    //корень 1 -
        {-1.0, -2.0},  //корень 2 -
        {-1.5, 1.5}    //корень 3 +
    };

    x1=-1.318459, x2=1.555637 
    x1=-0.881570, x2=-2.089052 
    x1=2.423889, x2=2.532399 
    x1=2.987711, x2=-0.916494 
    (0.350, −0.830)
    
    double solutions[][2] = {
        {2.0, 2.5},    //корень 1 - 
        {-1.0, -2.0},  //корень 2 -
        {-1.5, 1.5},    //корень 3 +
        {1.3, -1.8}
    };
    */
    
    
    
    double solutions[][2] = {
        {-1.5, 1.5},    // Для корня 1 (-1.318459, 1.555637)
        {-1.0, -2.0},  //корень 2 -
        {-1.5, 1.5},    //корень 3 +
        {3.0, -1.0}     // Для корня 4 (2.987711, -0.916494)
    };

    for (int i = 0; i < 4; i++) { //3
        printf("\n=== Корень %d ===\n", i+1);
        
        double x1_n = solutions[i][0], x2_n = solutions[i][1];
        double x1_s = solutions[i][0], x2_s = solutions[i][1];
        double x1_z = solutions[i][0], x2_z = solutions[i][1];
        
        printf("\n1. ");
        int iter_n = newton_method(&x1_n, &x2_n);
        
        printf("\n2. ");
        int iter_s = simple_iteration_adaptive(&x1_s, &x2_s);
        
        printf("\n3. ");
        int iter_z = seidel_method_adaptive(&x1_z, &x2_z);
        
        printf("\nИтоговые решения:\n");
        printf("Ньютон:    x1=%.6f, x2=%.6f (итераций: %d)\n", x1_n, x2_n, iter_n);
        if (iter_s > 0) printf("П.итерации: x1=%.6f, x2=%.6f (итераций: %d)\n", x1_s, x2_s, iter_s);
        if (iter_z > 0) printf("Зейделя:    x1=%.6f, x2=%.6f (итераций: %d)\n", x1_z, x2_z, iter_z);
    }
    
    return 0;
}