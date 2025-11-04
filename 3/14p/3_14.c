#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double f(double x) {
    return (log(x)*log(x) + sin(x/2.0)) / sqrt(x);
}

double df(double x) {
    // (u/v)' = (u'v - uv')/v²
    double u = log(x)*log(x) + sin(x/2.0);
    double du = (2*log(x))/x + (cos(x/2.0))/2.0;
    double v = sqrt(x);
    double dv = 1.0/(2.0*sqrt(x));
    
    return (du*v - u*dv) / (v*v);
}

double d2f(double x) {
    double h = 1e-2; 
    double h2 = h * h;
    
    if (h2 < 1e-15) {
        printf("Слишком маленький h²\n");
        return 0.0;
    }
    
    double result = (f(x+h) - 2.0*f(x) + f(x-h)) / h2;
    
    if (!isfinite(result) || fabs(result) > 1e10) {
        printf("Переполнение в d2f(x=%.2f): %.2e\n", x, result);
        return 0.0;
    }
    
    return result;
}

double d4f(double x) {
    double h = 5e-2; 
    double h4 = h * h * h * h;
    
    if (h4 < 1e-15) {
        return 0.0;
    }
    
    double numerator = f(x+2*h) - 4*f(x+h) + 6*f(x) - 4*f(x-h) + f(x-2*h);
    double result = numerator / h4;
    
    if (!isfinite(result) || fabs(result) > 1e8) {
        printf("Переполнение в d4f(x=%.2f): %.2e\n", x, result);
        return 0.0;
    }
    
    return result;
}

double d5f(double x) {
    double h = 1e-1; 
    double h5 = h * h * h * h * h;
    
    if (h5 < 1e-15) {
        return 0.0;
    }
    
    double numerator = -f(x+3*h) + 12*f(x+2*h) - 39*f(x+h) + 56*f(x) 
                     - 39*f(x-h) + 12*f(x-2*h) - f(x-3*h);
    double result = numerator / (6 * h5);
    
    if (!isfinite(result) || fabs(result) > 1e6) {
        printf("Переполнение в d5f(x=%.2f): %.2e\n", x, result);
        return 0.0;
    }
    
    return result;
}

//метод средних прямоугольников
double rectangle_middle(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    
    for (int i = 0; i < n; i++) {
        double x_i = a + i * h;
        double x_i1 = a + (i + 1) * h;
        double x_mid = (x_i + x_i1) / 2.0;
        sum += f(x_mid);
    }
    
    return h * sum;
}

//метод трапеций
double trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }
    
    return h * sum;
}

//метод Симпсона
double simpson(double a, double b, int n) {
    if (n % 2 != 0) n++; //убеждаемся, что n четное
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    //нечетные узлы
    for (int i = 1; i < n; i += 2) {
        double x = a + i * h;
        sum += 4.0 * f(x);
    }
    
    //четные узлы
    for (int i = 2; i < n; i += 2) {
        double x = a + i * h;
        sum += 2.0 * f(x);
    }
    
    return h * sum / 3.0;
}

//метод Эйлера
double euler(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }
    
    double correction = (h*h / 12.0) * (df(a) - df(b));
    
    return h * sum + correction;
}

double find_max_derivative(double a, double b, double (*df_func)(double), int samples) {
    double max_val = 0.0;
    double step = (b - a) / samples;
    
    for (int i = 0; i <= samples; i++) {
        double x = a + i * step;
        double deriv = fabs(df_func(x));
        if (deriv > max_val) {
            max_val = deriv;
        }
    }
    return max_val;
}

double error_rectangle_middle(double a, double b, int n) {
    double h = (b - a) / n;
    double max_f2 = find_max_derivative(a, b, d2f, 1000);
    return (b - a) * h * h * max_f2 / 24.0;
}

double error_trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double max_f2 = find_max_derivative(a, b, d2f, 1000);
    return (b - a) * h * h * max_f2 / 12.0;
}

double error_simpson(double a, double b, int n) {
    double h = (b - a) / n;
    double max_f4 = find_max_derivative(a, b, d4f, 1000);
    return (b - a) * h * h * h * h * max_f4 / 180.0;
}

double error_euler(double a, double b, int n) {
    double h = (b - a) / n;
    double max_f5 = find_max_derivative(a, b, d5f, 1000);
    return (b - a) * h * h * h * h * max_f5 / 720.0;
}

//уточнение по Рунге-Ромбергу
double runge_romberg(double I_h, double I_2h, int p) {
    return I_h + (I_h - I_2h) / (pow(2, p) - 1);
}

int main() {
    double a = 0.5, b = 2.0;
    int n1 = 8;  
    int n2 = 16; 
    
    printf("Вычисление интеграла: ∫[0.5, 2] (ln²x + sin(x/2)) / sqrt(x) dx\n\n");
    
    printf("=== Вычисление с шагом h (n = %d) ===\n", n1);
    
    double rect_h = rectangle_middle(a, b, n1);
    double trap_h = trapezoidal(a, b, n1);
    double simp_h = simpson(a, b, n1);
    double euler_h = euler(a, b, n1);
    
    printf("Метод средних прямоугольников: %.8f\n", rect_h);
    printf("Метод трапеций: %.8f\n", trap_h);
    printf("Метод Симпсона: %.8f\n", simp_h);
    printf("Метод Эйлера: %.8f\n\n", euler_h);
    
    printf("Оценки погрешностей (n = %d):\n", n1);
    printf("Средние прямоугольники: %.2f\n", error_rectangle_middle(a, b, n1));
    printf("Трапеций: %.2f\n", error_trapezoidal(a, b, n1));
    printf("Симпсона: %.2f\n", error_simpson(a, b, n1));
    printf("Эйлера: %.2f\n\n", error_euler(a, b, n1));
    
    printf("=== Вычисление с шагом h/2 (n = %d) ===\n", n2);
    
    double rect_2h = rectangle_middle(a, b, n2);
    double trap_2h = trapezoidal(a, b, n2);
    double simp_2h = simpson(a, b, n2);
    double euler_2h = euler(a, b, n2);
    
    printf("Метод средних прямоугольников: %.8f\n", rect_2h);
    printf("Метод трапеций: %.8f\n", trap_2h);
    printf("Метод Симпсона: %.8f\n", simp_2h);
    printf("Метод Эйлера: %.8f\n\n", euler_2h);
    
    printf("=== Уточнение по Рунге-Ромбергу ===\n");
    
    double rect_refined = runge_romberg(rect_2h, rect_h, 2);
    double trap_refined = runge_romberg(trap_2h, trap_h, 2);
    double simp_refined = runge_romberg(simp_2h, simp_h, 4);
    double euler_refined = runge_romberg(euler_2h, euler_h, 4);
    
    printf("Средние прямоугольники (уточненное): %.8f\n", rect_refined);
    printf("Трапеций (уточненное): %.8f\n", trap_refined);
    printf("Симпсона (уточненное): %.8f\n", simp_refined);
    printf("Эйлера (уточненное): %.8f\n\n", euler_refined);
    
    printf("=== Сводная таблица результатов ===\n");
    printf("Метод            | h (n=%d)   | h/2 (n=%d)| Уточненное\n", n1, n2);
    printf("--------------------------------------------------------\n");
    printf("Средние прямоуг. | %.6f  | %.6f  | %.6f\n", rect_h, rect_2h, rect_refined);
    printf("Трапеции         | %.6f  | %.6f  | %.6f\n", trap_h, trap_2h, trap_refined);
    printf("Симпсон          | %.6f  | %.6f  | %.6f\n", simp_h, simp_2h, simp_refined);
    printf("Эйлер            | %.6f  | %.6f  | %.6f\n", euler_h, euler_2h, euler_refined);
    
    return 0;
}