#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

// y = (5*ln(x^3 + 2) + x * e^x)/(2*x^2 + 3*x + 3)
double f(double x) {
    if (x <= 0) return 0;
    double numerator = 5 * log(x*x*x + 2) + x * exp(x);
    double denominator = 2*x*x + 3*x + 3;
    return numerator / denominator;
}

//аналитическая первая производная
double f_analytical_derivative(double x) {
    if (x <= 0) return 0;
    double u = 5 * log(x*x*x + 2) + x * exp(x);
    double v = 2*x*x + 3*x + 3;
    double u_prime = (15*x*x)/(x*x*x + 2) + exp(x) + x * exp(x);
    double v_prime = 4*x + 3;
    return (u_prime * v - u * v_prime) / (v * v);
}

//аналитическая вторая производная
double f_analytical_second_derivative(double x) {
    if (x <= 0) return 0;
    double u = 5 * log(x*x*x + 2) + x * exp(x);
    double v = 2*x*x + 3*x + 3;
    double u_prime = (15*x*x)/(x*x*x + 2) + exp(x) + x * exp(x);
    double v_prime = 4*x + 3;
    double u_double_prime = (30*x*(x*x*x + 2) - 15*x*x*3*x*x) / ((x*x*x + 2)*(x*x*x + 2)) + 2*exp(x) + x*exp(x);
    double v_double_prime = 4;
    double first_deriv = (u_prime * v - u * v_prime) / (v * v);
    double second_deriv = (u_double_prime*v - u*v_double_prime) / (v*v) - 2*(u_prime*v - u*v_prime)*v_prime/(v*v*v);
    return second_deriv;
}

// для первой производной
double right_difference(double x, double h, double (*func)(double)) {
    return (func(x + h) - func(x)) / h;
}

double left_difference(double x, double h, double (*func)(double)) {
    return (func(x) - func(x - h)) / h;
}

double central_difference(double x, double h, double (*func)(double)) {
    return (func(x + h) - func(x - h)) / (2 * h);
}

double three_point_forward(double x, double h, double (*func)(double)) {
    return (-3 * func(x) + 4 * func(x + h) - func(x + 2 * h)) / (2 * h);
}

double three_point_backward(double x, double h, double (*func)(double)) {
    return (3 * func(x) - 4 * func(x - h) + func(x - 2 * h)) / (2 * h);
}

double four_point(double x, double h, double (*func)(double)) {
    return (-2 * func(x - h) - 3 * func(x) + 6 * func(x + h) - func(x + 2 * h)) / (6 * h);
}

// для второй производной
double second_deriv_right(double x, double h, double (*func)(double)) {
    return (func(x) - 2 * func(x + h) + func(x + 2 * h)) / (h * h);
}

double second_deriv_left(double x, double h, double (*func)(double)) {
    return (func(x - 2 * h) - 2 * func(x - h) + func(x)) / (h * h);
}

double second_deriv_central(double x, double h, double (*func)(double)) {
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
}

double second_deriv_five_point(double x, double h, double (*func)(double)) {
    return (-func(x - 2 * h) + 16 * func(x - h) - 30 * func(x) + 16 * func(x + h) - func(x + 2 * h)) / (12 * h * h);
}

double absolute_error(double exact, double approx) {
    return fabs(exact - approx);
}

int main() {
    double a = 0.1;
    double b = 3.0;
    double h1 = 0.15;
    double h2 = 0.075;
    
    printf("ЧИСЛЕННОЕ ДИФФЕРЕНЦИРОВАНИЕ\n");
    printf("Функция: y = (5*ln(x^3 + 2) + x * e^x)/(2*x^2 + 3*x + 3)\n");
    printf("Интервал: [%.1f, %.1f]\n", a, b);
    printf("Шаги: h = %.3f и h/2 = %.3f\n\n", h1, h2);
    
    int n1 = (int)((b - a) / h1) + 1;
    int n2 = (int)((b - a) / h2) + 1;
    
    double *x1 = malloc(n1 * sizeof(double));
    double *x2 = malloc(n2 * sizeof(double));
    
    for (int i = 0; i < n1; i++) x1[i] = a + i * h1;
    for (int i = 0; i < n2; i++) x2[i] = a + i * h2;
    
    printf("ПЕРВАЯ ПРОИЗВОДНАЯ (h = %.3f)\n", h1);
    printf("============================================================================================================\n");
    printf("x\t\tAnalytic\tRight\t\tLeft\t\tCentral\t\t3pt_Fwd\t\t3pt_Bck\t\t4pt\n");
    printf("============================================================================================================\n");
    
    for (int i = 2; i < n1 - 2; i++) {
        double x = x1[i];
        double analytic = f_analytical_derivative(x);
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, analytic,
               right_difference(x, h1, f),
               left_difference(x, h1, f),
               central_difference(x, h1, f),
               three_point_forward(x, h1, f),
               three_point_backward(x, h1, f),
               four_point(x, h1, f));
    }
    
    printf("\nПЕРВАЯ ПРОИЗВОДНАЯ (h/2 = %.3f)\n", h2);
    printf("============================================================================================================\n");
    printf("x\t\tAnalytic\tRight\t\tLeft\t\tCentral\t\t3pt_Fwd\t\t3pt_Bck\t\t4pt\n");
    printf("============================================================================================================\n");
    
    for (int i = 4; i < n2 - 4; i += 2) {
        double x = x2[i];
        double analytic = f_analytical_derivative(x);
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, analytic,
               right_difference(x, h2, f),
               left_difference(x, h2, f),
               central_difference(x, h2, f),
               three_point_forward(x, h2, f),
               three_point_backward(x, h2, f),
               four_point(x, h2, f));
    }
    
    printf("\nВТОРАЯ ПРОИЗВОДНАЯ (h = %.3f)\n", h1);
    printf("======================================================================================\n");
    printf("x\t\tAnalytic\tRight\t\tLeft\t\tCentral\t\t5pt\n");
    printf("======================================================================================\n");
    
    for (int i = 2; i < n1 - 2; i++) {
        double x = x1[i];
        double analytic = f_analytical_second_derivative(x);
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, analytic,
               second_deriv_right(x, h1, f),
               second_deriv_left(x, h1, f),
               second_deriv_central(x, h1, f),
               second_deriv_five_point(x, h1, f));
    }
    
    printf("\nВТОРАЯ ПРОИЗВОДНАЯ (h/2 = %.3f)\n", h2);
    printf("======================================================================================\n");
    printf("x\t\tAnalytic\tRight\t\tLeft\t\tCentral\t\t5pt\n");
    printf("======================================================================================\n");
    
    for (int i = 4; i < n2 - 4; i += 2) {
        double x = x2[i];
        double analytic = f_analytical_second_derivative(x);
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, analytic,
               second_deriv_right(x, h2, f),
               second_deriv_left(x, h2, f),
               second_deriv_central(x, h2, f),
               second_deriv_five_point(x, h2, f));
    }
    
    printf("\nАНАЛИЗ ОШИБОК (h = %.3f)\n", h1);
    printf("============================================================================================================\n");
    printf("x\t\t1st_Err_Right\t1st_Err_Left\t1st_Err_Cent\t2nd_Err_Right\t2nd_Err_Left\t2nd_Err_Cent\n");
    printf("============================================================================================================\n");
    
    for (int i = 2; i < n1 - 2; i++) {
        double x = x1[i];
        double analytic1 = f_analytical_derivative(x);
        double analytic2 = f_analytical_second_derivative(x);
        
        double err1_right = absolute_error(analytic1, right_difference(x, h1, f));
        double err1_left = absolute_error(analytic1, left_difference(x, h1, f));
        double err1_cent = absolute_error(analytic1, central_difference(x, h1, f));
        
        double err2_right = absolute_error(analytic2, second_deriv_right(x, h1, f));
        double err2_left = absolute_error(analytic2, second_deriv_left(x, h1, f));
        double err2_cent = absolute_error(analytic2, second_deriv_central(x, h1, f));
        
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, err1_right, err1_left, err1_cent, err2_right, err2_left, err2_cent);
    }
    
    printf("\nАНАЛИЗ ОШИБОК (h/2 = %.3f)\n", h2);
    printf("============================================================================================================\n");
    printf("x\t\t1st_Err_Right\t1st_Err_Left\t1st_Err_Cent\t2nd_Err_Right\t2nd_Err_Left\t2nd_Err_Cent\n");
    printf("============================================================================================================\n");
    
    for (int i = 4; i < n2 - 4; i += 2) {
        double x = x2[i];
        double analytic1 = f_analytical_derivative(x);
        double analytic2 = f_analytical_second_derivative(x);
        
        double err1_right = absolute_error(analytic1, right_difference(x, h2, f));
        double err1_left = absolute_error(analytic1, left_difference(x, h2, f));
        double err1_cent = absolute_error(analytic1, central_difference(x, h2, f));
        
        double err2_right = absolute_error(analytic2, second_deriv_right(x, h2, f));
        double err2_left = absolute_error(analytic2, second_deriv_left(x, h2, f));
        double err2_cent = absolute_error(analytic2, second_deriv_central(x, h2, f));
        
        printf("%.3f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               x, err1_right, err1_left, err1_cent, err2_right, err2_left, err2_cent);
    }
    
    free(x1);
    free(x2);
    
    return 0;
}