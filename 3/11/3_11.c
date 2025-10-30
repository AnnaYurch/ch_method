#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void natural_cubic_spline(int n, double x[], double y[], double a[], double b[], double c[], double d[]) {
    double *h = (double*)malloc(n * sizeof(double)); // длина i-го отрезка между точками
    double *alpha = (double*)malloc(n * sizeof(double));
    double *l = (double*)malloc((n+1) * sizeof(double));
    double *mu = (double*)malloc((n+1) * sizeof(double));
    double *z = (double*)malloc((n+1) * sizeof(double));
    
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    c[0] = 0.0;
    c[n] = 0.0;
    
    //решаем трёхдиагональную систему для c[i]
    for (int i = 1; i < n; i++) {
        alpha[i] = 3.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
    }
    
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    
    //метод прогонки
    for (int i = 1; i < n; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }
    
    l[n] = 1.0;
    z[n] = 0.0;
    
    for (int j = n-1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j+1];
    }
    
    for (int i = 0; i < n; i++) {
        a[i] = y[i];
        b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (c[i+1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i]);
    }
    
    free(h);
    free(alpha);
    free(l);
    free(mu);
    free(z);
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

// Исправленная функция для построения графика
void plot_spline(int n, double x[], double y[], double a[], double b[], double c[], double d[], double x_star) {
    FILE *data_file = fopen("spline_data.txt", "w");
    FILE *gnuplot_script = fopen("plot_spline.gnu", "w");
    
    if (!data_file || !gnuplot_script) {
        printf("Ошибка создания файлов для графиков\n");
        return;
    }
    
    // Записываем исходные точки
    fprintf(data_file, "# ИСХОДНЫЕ ТОЧКИ\n");
    for (int i = 0; i <= n; i++) {
        fprintf(data_file, "%.6f %.6f\n", x[i], y[i]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точки сплайна БЕЗ разрывов между отрезками
    fprintf(data_file, "# КУБИЧЕСКИЙ СПЛАЙН\n");
    for (int i = 0; i < n; i++) {
        // Вычисляем точки на каждом отрезке
        int points_per_segment = 20;
        for (int j = 0; j <= points_per_segment; j++) {
            double xi = x[i] + j * (x[i+1] - x[i]) / points_per_segment;
            double yi = spline_value(n, x, a, b, c, d, xi);
            fprintf(data_file, "%.6f %.6f\n", xi, yi);
        }
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точку интерполяции
    fprintf(data_file, "# ТОЧКА ИНТЕРПОЛЯЦИИ\n");
    double y_star = spline_value(n, x, a, b, c, d, x_star);
    fprintf(data_file, "%.6f %.6f\n", x_star, y_star);
    
    fclose(data_file);
    
    // Вычисляем min и max для установки диапазона Y
    double y_min = min_element(y, n+1);
    double y_max = max_element(y, n+1);
    
    // Создаем скрипт для gnuplot
    fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
    fprintf(gnuplot_script, "set output 'spline_graph.png'\n");
    fprintf(gnuplot_script, "set title 'Естественный кубический сплайн дефекта 1' font 'Arial,14'\n");
    fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set ylabel 'y' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set grid\n");
    fprintf(gnuplot_script, "set key top left box\n");
    fprintf(gnuplot_script, "set xrange [%.1f:%.1f]\n", x[0]-0.1, x[n]+0.1);
    fprintf(gnuplot_script, "set yrange [%.1f:%.1f]\n", y_min - 0.5, y_max + 0.5);
    
    // Настройка стилей
    fprintf(gnuplot_script, "set style line 1 lc rgb 'black' pt 7 ps 1.5\n");
    fprintf(gnuplot_script, "set style line 2 lc rgb 'red' lw 2\n");
    fprintf(gnuplot_script, "set style line 3 lc rgb 'blue' pt 9 ps 2\n");
    
    fprintf(gnuplot_script, "plot 'spline_data.txt' index 0 with points ls 1 title 'Исходные точки', \\\n");
    fprintf(gnuplot_script, "     'spline_data.txt' index 1 with lines ls 2 title 'Кубический сплайн', \\\n");
    fprintf(gnuplot_script, "     'spline_data.txt' index 2 with points ls 3 title 'x*=%.3f (S=%.3f)'\n", 
            x_star, y_star);
    
    fclose(gnuplot_script);
    
    // Запускаем gnuplot
    system("gnuplot plot_spline.gnu");
    printf("График сохранен в файл: spline_graph.png\n");
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
    printf("\n");
    
    printf("ВСЕ КОЭФФИЦИЕНТЫ СПЛАЙНА:\n");
    printf("i\tx[i]\ty[i]\ta[i]\tb[i]\tc[i]\td[i]\n");
    for (int i = 0; i < n; i++) {
        printf("%d\t%.3f\t%.3f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
               i, x[i], y[i], a[i], b[i], c[i], d[i]);
    }
    printf("%d\t%.3f\t%.3f\t\t\t%.6f\n", n, x[n], y[n], c[n]);
    
    printf("\n=== ПОСТРОЕНИЕ ГРАФИКА ===\n");
    plot_spline(n, x, y, a, b, c, d, x_star);

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}