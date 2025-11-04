#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Функция для вычисления ПЕРВОЙ ПРОИЗВОДНОЙ сплайна в точке xp
double spline_derivative(int n, double x[], double a[], double b[], double c[], double d[], double xp) {
    for (int i = 0; i < n; i++) {
        if (xp >= x[i] && xp <= x[i+1]) {
            double dx = xp - x[i];
            // Производная полинома: S'(x) = b_i + 2*c_i*dx + 3*d_i*dx^2
            return b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;
        }
    }
    return 0.0;
}

// Функция для вычисления ВТОРОЙ ПРОИЗВОДНОЙ сплайна в точке xp
double spline_second_derivative(int n, double x[], double a[], double b[], double c[], double d[], double xp) {
    for (int i = 0; i < n; i++) {
        if (xp >= x[i] && xp <= x[i+1]) {
            double dx = xp - x[i];
            // Вторая производная полинома: S''(x) = 2*c_i + 6*d_i*dx
            return 2.0 * c[i] + 6.0 * d[i] * dx;
        }
    }
    return 0.0;
}

// Классический метод прогонки по формулам из PDF
void classical_sweep_method(int n, double a[], double b[], double c[], double d_vec[], double x[]) {
    double *P = (double*)malloc(n * sizeof(double));
    double *Q = (double*)malloc(n * sizeof(double));
    
    // ПРЯМОЙ ХОД (вычисление прогоночных коэффициентов)
    
    // Первое уравнение: b₀x₀ + c₀x₁ = d₀
    P[0] = -c[0] / b[0];
    Q[0] = d_vec[0] / b[0];
    
    // Промежуточные уравнения: aᵢxᵢ₋₁ + bᵢxᵢ + cᵢxᵢ₊₁ = dᵢ
    for (int i = 1; i < n-1; i++) {
        double denominator = b[i] + a[i] * P[i-1];
        P[i] = -c[i] / denominator;
        Q[i] = (d_vec[i] - a[i] * Q[i-1]) / denominator;
    }
    
    // Последнее уравнение: aₙ₋₁xₙ₋₂ + bₙ₋₁xₙ₋₁ = dₙ₋₁
    int last = n-1;
    double denominator = b[last] + a[last] * P[last-1];
    P[last] = 0.0;  // По формуле: Pₙ = 0
    Q[last] = (d_vec[last] - a[last] * Q[last-1]) / denominator;
    
    // ОБРАТНЫЙ ХОД (вычисление решения)
    x[last] = Q[last];
    for (int i = n-2; i >= 0; i--) {
        x[i] = P[i] * x[i+1] + Q[i];
    }
    
    free(P);
    free(Q);
}

// Построение естественного кубического сплайна с использованием классического метода прогонки
void natural_cubic_spline(int n, double x[], double y[], double a[], double b[], double c[], double d[]) {
    double *h = (double*)malloc(n * sizeof(double));
    
    // Вычисляем длины отрезков
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    // Естественные краевые условия
    c[0] = 0.0;
    c[n] = 0.0;
    
    // Если у нас только 2 точки (n=1), то сплайн - прямая линия
    if (n == 1) {
        a[0] = y[0];
        b[0] = (y[1] - y[0]) / h[0];
        c[0] = 0.0;
        d[0] = 0.0;
        free(h);
        return;
    }
    
    // Формируем систему для c[1]...c[n-1]
    int system_size = n - 1; // количество уравнений в системе
    
    double *A = (double*)malloc(system_size * sizeof(double)); // нижняя диагональ (a)
    double *B = (double*)malloc(system_size * sizeof(double)); // главная диагональ (b)  
    double *C = (double*)malloc(system_size * sizeof(double)); // верхняя диагональ (c)
    double *D = (double*)malloc(system_size * sizeof(double)); // правая часть (d)
    double *c_solution = (double*)malloc(system_size * sizeof(double)); // решения (c_i)
    
    // Заполняем коэффициенты системы
    for (int i = 0; i < system_size; i++) {
        int node_index = i + 1; // соответствие: i в системе -> узел i+1 в сплайне
        
        if (i == 0) {
            // Первое уравнение системы (для узла 1)
            A[i] = 0; // a₀ = 0
            B[i] = 2.0 * (h[0] + h[1]); // b₀
            C[i] = h[1]; // c₀
            D[i] = 3.0 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0]);
        }
        else if (i == system_size - 1) {
            // Последнее уравнение системы (для узла n-1)
            A[i] = h[n-2]; // a_{n-2}
            B[i] = 2.0 * (h[n-2] + h[n-1]); // b_{n-2}
            C[i] = 0; // c_{n-2} = 0
            D[i] = 3.0 * ((y[n] - y[n-1]) / h[n-1] - (y[n-1] - y[n-2]) / h[n-2]);
        }
        else {
            // Промежуточные уравнения
            A[i] = h[node_index - 1]; // aᵢ
            B[i] = 2.0 * (h[node_index - 1] + h[node_index]); // bᵢ
            C[i] = h[node_index]; // cᵢ
            D[i] = 3.0 * ((y[node_index + 1] - y[node_index]) / h[node_index] - 
                          (y[node_index] - y[node_index - 1]) / h[node_index - 1]);
        }
    }
    
    // Решаем систему методом прогонки
    classical_sweep_method(system_size, A, B, C, D, c_solution);
    
    // Переносим решения в массив c[]
    for (int i = 0; i < system_size; i++) {
        c[i + 1] = c_solution[i];
    }
    
    // Вычисляем остальные коэффициенты сплайна
    for (int i = 0; i < n; i++) {
        a[i] = y[i];
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
    }
    
    free(h);
    free(A);
    free(B);
    free(C);
    free(D);
    free(c_solution);
}

double spline_value(int n, double x[], double a[], double b[], double c[], double d[], double xp) {
    for (int i = 0; i < n; i++) {
        if (xp >= x[i] && xp <= x[i+1]) {
            double dx = xp - x[i];
            return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    }
    return 0.0; 
}

double min_element(double arr[], int size) {
    double min_val = arr[0];
    for (int i = 1; i < size; i++) {
        if (arr[i] < min_val) min_val = arr[i];
    }
    return min_val;
}

double max_element(double arr[], int size) {
    double max_val = arr[0];
    for (int i = 1; i < size; i++) {
        if (arr[i] > max_val) max_val = arr[i];
    }
    return max_val;
}

// Функция для проверки непрерывности производных в узлах
// Функция для проверки непрерывности производных в узлах
void check_continuity(int n, double x[], double a[], double b[], double c[], double d[]) {
    printf("\n=== ПРОВЕРКА НЕПРЕРЫВНОСТИ ПРОИЗВОДНЫХ ===\n");
    
    // Сначала вычисляем длины отрезков для проверки второй производной
    double *h = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    // ПРОВЕРКА ВТОРОЙ ПРОИЗВОДНОЙ
    printf("Проверка второй производной в узлах:\n");
    for (int i = 1; i < n; i++) {
        double left = 2*c[i-1] + 6*d[i-1]*h[i-1];  // S'' слева
        double right = 2*c[i];                     // S'' справа  
        printf("Узел %d (x=%.3f): S''(слева)=%.6f, S''(справа)=%.6f, разность=%.6f\n", 
               i, x[i], left, right, fabs(left-right));
    }
    printf("\n");
    
    // ПРОВЕРКА ПЕРВОЙ ПРОИЗВОДНОЙ (существующий код)
    for (int i = 1; i < n; i++) {
        // Левая производная в узле x[i] (от отрезка i-1)
        double dx_left = x[i] - x[i-1];
        double S_left = b[i-1] + 2.0 * c[i-1] * dx_left + 3.0 * d[i-1] * dx_left * dx_left;
        
        // Правая производная в узле x[i] (от отрезка i)
        double S_right = b[i];
        
        // Вторая производная слева
        double S2_left = 2.0 * c[i-1] + 6.0 * d[i-1] * dx_left;
        
        // Вторая производная справа
        double S2_right = 2.0 * c[i];
        
        printf("Узел x[%d] = %.3f:\n", i, x[i]);
        printf("  S'(слева) = %.6f, S'(справа) = %.6f, разность = %.6f\n", 
               S_left, S_right, fabs(S_left - S_right));
        printf("  S''(слева) = %.6f, S''(справа) = %.6f, разность = %.6f\n", 
               S2_left, S2_right, fabs(S2_left - S2_right));
    }
    
    free(h);
}

// Функция для построения графика ТОЛЬКО первой производной
void plot_first_derivative(int n, double x[], double y[], double a[], double b[], double c[], double d[], double x_star) {
    FILE *data_file = fopen("derivative_data.txt", "w");
    
    if (!data_file) {
        printf("Ошибка создания файла данных для графиков\n");
        return;
    }
    
    // Записываем исходные точки (для контекста)
    fprintf(data_file, "# ИСХОДНЫЕ ТОЧКИ\n");
    for (int i = 0; i <= n; i++) {
        fprintf(data_file, "%.6f %.6f\n", x[i], y[i]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точки ПЕРВОЙ ПРОИЗВОДНОЙ сплайна
    fprintf(data_file, "# ПЕРВАЯ ПРОИЗВОДНАЯ СПЛАЙНА\n");
    for (int i = 0; i < n; i++) {
        int points_per_segment = 50; // Увеличиваем количество точек для гладкого графика
        for (int j = 0; j <= points_per_segment; j++) {
            double xi = x[i] + j * (x[i+1] - x[i]) / points_per_segment;
            double dyi = spline_derivative(n, x, a, b, c, d, xi);
            fprintf(data_file, "%.6f %.6f\n", xi, dyi);
        }
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точку интерполяции
    fprintf(data_file, "# ТОЧКА ИНТЕРПОЛЯЦИИ\n");
    double y_derivative = spline_derivative(n, x, a, b, c, d, x_star);
    fprintf(data_file, "%.6f %.6f\n", x_star, y_derivative);
    
    fclose(data_file);
    
    // Вычисляем диапазон для первой производной
    double dy_min = 1e10, dy_max = -1e10;
    for (int i = 0; i < n; i++) {
        int points_per_segment = 10;
        for (int j = 0; j <= points_per_segment; j++) {
            double xi = x[i] + j * (x[i+1] - x[i]) / points_per_segment;
            double dyi = spline_derivative(n, x, a, b, c, d, xi);
            if (dyi < dy_min) dy_min = dyi;
            if (dyi > dy_max) dy_max = dyi;
        }
    }
    
    // Добавляем немного запаса по краям
    dy_min -= 0.5;
    dy_max += 0.5;
    
    // Создаем скрипт для gnuplot (ТОЛЬКО первая производная)
    FILE *gnuplot_script = fopen("plot_first_derivative.gnu", "w");
    if (gnuplot_script) {
        fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
        fprintf(gnuplot_script, "set output 'first_derivative.png'\n");
        fprintf(gnuplot_script, "set title 'Первая производная кубического сплайна' font 'Arial,14'\n");
        fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
        fprintf(gnuplot_script, "set ylabel 'S''(x)' font 'Arial,12'\n");
        fprintf(gnuplot_script, "set grid\n");
        fprintf(gnuplot_script, "set key top left box\n");
        fprintf(gnuplot_script, "set xrange [%.1f:%.1f]\n", x[0]-0.1, x[n]+0.1);
        fprintf(gnuplot_script, "set yrange [%.1f:%.1f]\n", dy_min, dy_max);
        
        // Настройка стилей
        fprintf(gnuplot_script, "set style line 1 lc rgb 'blue' lw 3\n");
        fprintf(gnuplot_script, "set style line 2 lc rgb 'red' pt 9 ps 2\n");
        fprintf(gnuplot_script, "set style line 3 lc rgb 'black' pt 7 ps 1.5\n");
        
        fprintf(gnuplot_script, "plot 'derivative_data.txt' index 1 with lines ls 1 title 'Первая производная S'(x)', \\\n");
        fprintf(gnuplot_script, "     'derivative_data.txt' index 2 with points ls 2 title 'x*=%.3f (S'=%.3f)', \\\n", 
                x_star, y_derivative);
        fprintf(gnuplot_script, "     'derivative_data.txt' index 0 with points ls 3 title 'Исходные точки'\n");
        
        fclose(gnuplot_script);
        
        // Запускаем gnuplot
        system("gnuplot plot_first_derivative.gnu");
        printf("График первой производной сохранен в файл: first_derivative.png\n");
    }
}

// Функция для построения графика ТОЛЬКО второй производной
void plot_second_derivative(int n, double x[], double y[], double a[], double b[], double c[], double d[], double x_star) {
    FILE *data_file = fopen("second_derivative_data.txt", "w");
    
    if (!data_file) {
        printf("Ошибка создания файла данных для графиков\n");
        return;
    }
    
    // Записываем исходные точки (для контекста)
    fprintf(data_file, "# ИСХОДНЫЕ ТОЧКИ\n");
    for (int i = 0; i <= n; i++) {
        fprintf(data_file, "%.6f %.6f\n", x[i], y[i]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точки ВТОРОЙ ПРОИЗВОДНОЙ сплайна
    fprintf(data_file, "# ВТОРАЯ ПРОИЗВОДНАЯ СПЛАЙНА\n");
    for (int i = 0; i < n; i++) {
        int points_per_segment = 100; // Увеличиваем количество точек для гладкого графика
        for (int j = 0; j <= points_per_segment; j++) {
            double xi = x[i] + j * (x[i+1] - x[i]) / points_per_segment;
            double d2yi = spline_second_derivative(n, x, a, b, c, d, xi);
            fprintf(data_file, "%.6f %.6f\n", xi, d2yi);
        }
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точку интерполяции
    fprintf(data_file, "# ТОЧКА ИНТЕРПОЛЯЦИИ\n");
    double y_second_derivative = spline_second_derivative(n, x, a, b, c, d, x_star);
    fprintf(data_file, "%.6f %.6f\n", x_star, y_second_derivative);
    
    fclose(data_file);
    
    // Вычисляем диапазон для второй производной
    double d2y_min = 1e10, d2y_max = -1e10;
    for (int i = 0; i < n; i++) {
        int points_per_segment = 20;
        for (int j = 0; j <= points_per_segment; j++) {
            double xi = x[i] + j * (x[i+1] - x[i]) / points_per_segment;
            double d2yi = spline_second_derivative(n, x, a, b, c, d, xi);
            if (d2yi < d2y_min) d2y_min = d2yi;
            if (d2yi > d2y_max) d2y_max = d2yi;
        }
    }
    
    // Добавляем запас по краям и убеждаемся, что 0 включен в диапазон
    double range = d2y_max - d2y_min;
    d2y_min = fmin(d2y_min - 0.1 * range, -5.0);
    d2y_max = fmax(d2y_max + 0.1 * range, 5.0);
    
    // Создаем скрипт для gnuplot (ТОЛЬКО вторая производная)
    FILE *gnuplot_script = fopen("plot_second_derivative.gnu", "w");
    if (gnuplot_script) {
        fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
        fprintf(gnuplot_script, "set output 'second_derivative.png'\n");
        fprintf(gnuplot_script, "set title 'Вторая производная кубического сплайна' font 'Arial,14'\n");
        fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
        fprintf(gnuplot_script, "set ylabel 'S''''(x)' font 'Arial,12'\n");
        fprintf(gnuplot_script, "set grid\n");
        fprintf(gnuplot_script, "set key top left box\n");
        fprintf(gnuplot_script, "set xrange [%.1f:%.1f]\n", x[0]-0.1, x[n]+0.1);
        fprintf(gnuplot_script, "set yrange [%.1f:%.1f]\n", d2y_min, d2y_max);
        
        // Настройка стилей
        fprintf(gnuplot_script, "set style line 1 lc rgb 'purple' lw 3\n");
        fprintf(gnuplot_script, "set style line 2 lc rgb 'orange' pt 9 ps 2\n");
        fprintf(gnuplot_script, "set style line 3 lc rgb 'black' pt 7 ps 1.5\n");
        
        fprintf(gnuplot_script, "plot 'second_derivative_data.txt' index 1 with lines ls 1 title 'Вторая производная', \\\n");
        fprintf(gnuplot_script, "     'second_derivative_data.txt' index 2 with points ls 2 title 'x*=%.3f (d²S/dx²=%.3f)', \\\n", 
                x_star, y_second_derivative);
        fprintf(gnuplot_script, "     'second_derivative_data.txt' index 0 with points ls 3 title 'Исходные точки'\n");
        
        fclose(gnuplot_script);
        
        // Запускаем gnuplot
        system("gnuplot plot_second_derivative.gnu");
        printf("График второй производной сохранен в файл: second_derivative.png\n");
    }
}

int main() {
    int n = 10; //количество отрезков
    double x[] = {0.06, 0.356, 0.763, 1.096, 1.54, 2.021, 2.354, 2.761, 3.205, 3.464, 3.76};
    double y[] = {2.478, 3.471, 3.257, 1.086, 1.291, 3.735, 0.729, 0.684, 2.371, 3.142, 2.846};
    
    double *a = (double*)malloc(n * sizeof(double));
    double *b = (double*)malloc(n * sizeof(double));
    double *c = (double*)malloc((n+1) * sizeof(double));
    double *d = (double*)malloc(n * sizeof(double));
    
    natural_cubic_spline(n, x, y, a, b, c, d);
    
    double x_star = 1.285;
    
    //находим отрезок для x_star
    int segment = -1;
    for (int i = 0; i < n; i++) {
        if (x_star >= x[i] && x_star <= x[i+1]) {
            segment = i;
            break;
        }
    }
    
    if (segment == -1) {
        printf("Точка x* = %.3f вне диапазона данных\n", x_star);
        return 1;
    }
    
    //вычисляем значение в точке x_star
    double y_star = spline_value(n, x, a, b, c, d, x_star);
    double y_derivative = spline_derivative(n, x, a, b, c, d, x_star);
    double y_second_derivative = spline_second_derivative(n, x, a, b, c, d, x_star);
    
    printf("=== ЕСТЕСТВЕННЫЙ КУБИЧЕСКИЙ СПЛАЙН ДЕФЕКТА 1 ===\n");
    printf("Точка x* = %.3f находится на отрезке [%.3f, %.3f] (i = %d)\n", 
           x_star, x[segment], x[segment+1], segment+1);
    printf("Коэффициенты сплайна на этом отрезке:\n");
    printf("a[%d] = %.6f\n", segment+1, a[segment]);
    printf("b[%d] = %.6f\n", segment+1, b[segment]);
    printf("c[%d] = %.6f\n", segment+1, c[segment]);
    printf("d[%d] = %.6f\n", segment+1, d[segment]);
    printf("\n");
    printf("Значение сплайна в точке x* = %.3f:\n", x_star);
    printf("S(%.3f) = %.6f\n", x_star, y_star);
    printf("S'(%.3f) = %.6f\n", x_star, y_derivative);
    printf("S''(%.3f) = %.6f\n", x_star, y_second_derivative);
    printf("\n");
    
    printf("ВСЕ КОЭФФИЦИЕНТЫ СПЛАЙНА:\n");
    printf("i\tx[i]\ty[i]\ta[i]\tb[i]\tc[i]\td[i]\n");
    for (int i = 0; i < n; i++) {
        printf("%d\t%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               i, x[i], y[i], a[i], b[i], c[i], d[i]);
    }
    printf("%d\t%.3f\t%.3f\t\t\t%.6f\n", n, x[n], y[n], c[n]);
    
    // Проверка непрерывности производных
    check_continuity(n, x, a, b, c, d);
    
    printf("\n=== ПОСТРОЕНИЕ ГРАФИКА ПЕРВОЙ ПРОИЗВОДНОЙ ===\n");
    //plot_first_derivative(n, x, y, a, b, c, d, x_star);
    plot_second_derivative(n, x, y, a, b, c, d, x_star);

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}