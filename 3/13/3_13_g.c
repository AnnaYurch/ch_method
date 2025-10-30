#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

// Простая функция для создания графика
void create_simple_first_derivative_plot(double a, double b, double h) {
    printf("Создание данных для первой производной (h=%.3f)...\n", h);
    
    FILE *data_file = fopen("first_deriv_simple.txt", "w");
    if (!data_file) {
        printf("Ошибка создания файла данных\n");
        return;
    }
    
    // Сужаем интервал для избежания граничных эффектов
    double plot_a = a + 2*h;
    double plot_b = b - 2*h;
    
    int n_points = 30;
    double step = (plot_b - plot_a) / n_points;
    
    // Записываем заголовок
    fprintf(data_file, "# x Analytic Right Left Central\n");
    
    // Записываем данные только для основных методов
    for (int i = 0; i <= n_points; i++) {
        double x = plot_a + i * step;
        double analytic = f_analytical_derivative(x);
        double right = right_difference(x, h, f);
        double left = left_difference(x, h, f);
        double central = central_difference(x, h, f);
        
        fprintf(data_file, "%.6f %.6f %.6f %.6f %.6f\n", 
               x, analytic, right, left, central);
    }
    fclose(data_file);
    
    printf("Данные записаны в first_deriv_simple.txt\n");
    
    // Создаем простой скрипт gnuplot
    FILE *gnuplot_script = fopen("plot_first_simple.gnu", "w");
    if (!gnuplot_script) {
        printf("Ошибка создания скрипта gnuplot\n");
        return;
    }
    
    fprintf(gnuplot_script, "set terminal pngcairo size 800,600\n");
    fprintf(gnuplot_script, "set output 'first_derivative_simple.png'\n");
    fprintf(gnuplot_script, "set title 'First Derivative Methods (h=%.3f)'\n", h);
    fprintf(gnuplot_script, "set xlabel 'x'\n");
    fprintf(gnuplot_script, "set ylabel 'f'(x)'\n");
    fprintf(gnuplot_script, "set grid\n");
    fprintf(gnuplot_script, "set key top left\n");
    
    fprintf(gnuplot_script, "plot 'first_deriv_simple.txt' using 1:2 with lines lw 3 title 'Analytic', \\\n");
    fprintf(gnuplot_script, "     '' using 1:3 with lines lw 2 title 'Right Difference', \\\n");
    fprintf(gnuplot_script, "     '' using 1:4 with lines lw 2 title 'Left Difference', \\\n");
    fprintf(gnuplot_script, "     '' using 1:5 with lines lw 2 title 'Central Difference'\n");
    
    fclose(gnuplot_script);
    
    printf("Запуск gnuplot...\n");
    int result = system("gnuplot plot_first_simple.gnu");
    
    if (result == 0) {
        printf("УСПЕХ: График создан - first_derivative_simple.png\n");
    } else {
        printf("ОШИБКА: Не удалось создать график. Проверьте установку gnuplot.\n");
    }
}

// Функция для проверки данных
void check_data_and_system() {
    printf("\n=== ПРОВЕРКА СИСТЕМЫ ===\n");
    
    // Проверяем gnuplot
    printf("1. Проверка gnuplot: ");
    int gnuplot_check = system("which gnuplot > /dev/null 2>&1");
    if (gnuplot_check == 0) {
        printf("УСТАНОВЛЕН\n");
    } else {
        printf("НЕ УСТАНОВЛЕН\n");
    }
    
    // Проверяем создание простого тестового графика
    printf("2. Создание тестового графика...\n");
    FILE *test_script = fopen("test_plot.gnu", "w");
    if (test_script) {
        fprintf(test_script, "set terminal png\n");
        fprintf(test_script, "set output 'test_output.png'\n");
        fprintf(test_script, "plot sin(x)\n");
        fclose(test_script);
        
        system("gnuplot test_plot.gnu 2>/dev/null");
        
        // Проверяем создался ли файл
        FILE *test_file = fopen("test_output.png", "r");
        if (test_file) {
            printf("   ТЕСТ УСПЕШЕН: PNG файл создан\n");
            fclose(test_file);
        } else {
            printf("   ТЕСТ НЕ УДАЛСЯ: PNG файл не создан\n");
        }
    }
}

// Функция для текстового вывода графика
void print_text_plot(double a, double b, double h) {
    printf("\n=== ТЕКСТОВЫЙ ГРАФИК ПЕРВОЙ ПРОИЗВОДНОЙ (h=%.3f) ===\n", h);
    printf("x\tAnalytic\tRight\t\tLeft\t\tCentral\t\t3pt_Fwd\t\t3pt_Bck\t\t4pt\n");
    printf("--------------------------------------------------------------------------------------------------------\n");
    
    double plot_a = a + 2*h;
    double plot_b = b - 2*h;
    int n_points = 15;
    double step = (plot_b - plot_a) / n_points;
    
    for (int i = 0; i <= n_points; i++) {
        double x = plot_a + i * step;
        double analytic = f_analytical_derivative(x);
        double right = right_difference(x, h, f);
        double left = left_difference(x, h, f);
        double central = central_difference(x, h, f);
        double three_fwd = three_point_forward(x, h, f);
        double three_bck = three_point_backward(x, h, f);
        double four_pt = four_point(x, h, f);
        
        printf("%.2f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n", 
               x, analytic, right, left, central, three_fwd, three_bck, four_pt);
    }
}

int main() {
    double a = 0.1;
    double b = 3.0;
    double h1 = 0.15;
    double h2 = 0.075;
    
    printf("ЧИСЛЕННОЕ ДИФФЕРЕНЦИРОВАНИЕ - ПЕРВАЯ ПРОИЗВОДНАЯ\n");
    printf("=================================================\n");
    
    // Сначала проверяем систему
    check_data_and_system();
    
    // Выводим текстовый график
    printf("\n=== ТЕКСТОВОЕ ПРЕДСТАВЛЕНИЕ ===\n");
    print_text_plot(a, b, h1);
    
    // Пробуем создать простой график
    printf("\n=== СОЗДАНИЕ ГРАФИКА ===\n");
    create_simple_first_derivative_plot(a, b, h1);
    
    // Проверяем созданные файлы
    printf("\n=== ПРОВЕРКА СОЗДАННЫХ ФАЙЛОВ ===\n");
    system("ls -la first_deriv_simple.txt 2>/dev/null || echo 'Файл данных не найден'");
    system("ls -la first_derivative_simple.png 2>/dev/null || echo 'PNG файл не найден'");
    
    // Показываем первые строки данных
    printf("\n=== ПЕРВЫЕ 5 СТРОК ДАННЫХ ===\n");
    system("head -5 first_deriv_simple.txt 2>/dev/null || echo 'Файл данных не существует'");
    
    return 0;
}