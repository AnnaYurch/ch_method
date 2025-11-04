#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double lagrange(double x, double x_nodes[], double y_nodes[], int n, int indices[]) {
    double result = 0.0;
    
    for (int i = 0; i <= n; i++) {
        double term = y_nodes[indices[i]];
        for (int j = 0; j <= n; j++) {
            if (j != i) {
                term *= (x - x_nodes[indices[j]]) / 
                       (x_nodes[indices[i]] - x_nodes[indices[j]]);
            }
        }
        result += term;
    }
    
    return result;
}

//функция для проверки в узловой точке
void check_at_node(double x_nodes[], double y_nodes[], int node_index, int indices[], int n) {
    double interpolated = lagrange(x_nodes[node_index], x_nodes, y_nodes, n, indices);
    //значение многочлена Лагранжа в узловой точке
    double actual = y_nodes[node_index];
    double error = fabs(interpolated - actual);
    
    printf("Проверка в узле x=%.2f:\n", x_nodes[node_index]);
    printf("  Истинное значение: %.4f\n", actual);
    printf("  Интерполированное: %.4f\n", interpolated);
    printf("  Ошибка: %.6f\n\n", error);
}

void plot_lagrange(double x_nodes[], double y_nodes[], int indices_L2[], int indices_L3[], int n_total, double x_star) {
    FILE *data_file = fopen("lagrange_data.txt", "w");
    FILE *gnuplot_script = fopen("plot_lagrange.gnu", "w");
    
    if (!data_file || !gnuplot_script) {
        printf("Ошибка создания файлов для графиков\n");
        return;
    }
    
    // Записываем исходные точки
    fprintf(data_file, "# ИСХОДНЫЕ ТОЧКИ\n");
    for (int i = 0; i < n_total; i++) {
        fprintf(data_file, "%.6f %.6f\n", x_nodes[i], y_nodes[i]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точки для L2 (многочлен 2-й степени)
    fprintf(data_file, "# МНОГОЧЛЕН L2(x) 2-й СТЕПЕНИ\n");
    for (double xi = -1.5; xi <= 2.5; xi += 0.01) {
        double yi = lagrange(xi, x_nodes, y_nodes, 2, indices_L2);
        fprintf(data_file, "%.6f %.6f\n", xi, yi);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точки для L3 (многочлен 3-й степени)
    fprintf(data_file, "# МНОГОЧЛЕН L3(x) 3-й СТЕПЕНИ\n");
    for (double xi = -1.5; xi <= 2.5; xi += 0.01) {
        double yi = lagrange(xi, x_nodes, y_nodes, 3, indices_L3);
        fprintf(data_file, "%.6f %.6f\n", xi, yi);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем узлы, используемые для L2
    fprintf(data_file, "# УЗЛЫ L2\n");
    for (int i = 0; i < 3; i++) {
        fprintf(data_file, "%.6f %.6f\n", x_nodes[indices_L2[i]], y_nodes[indices_L2[i]]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем узлы, используемые для L3
    fprintf(data_file, "# УЗЛЫ L3\n");
    for (int i = 0; i < 4; i++) {
        fprintf(data_file, "%.6f %.6f\n", x_nodes[indices_L3[i]], y_nodes[indices_L3[i]]);
    }
    fprintf(data_file, "\n\n");
    
    // Записываем точку x_star
    fprintf(data_file, "# ТОЧКА X_STAR\n");
    double y_star_L2 = lagrange(x_star, x_nodes, y_nodes, 2, indices_L2);
    double y_star_L3 = lagrange(x_star, x_nodes, y_nodes, 3, indices_L3);
    fprintf(data_file, "%.6f %.6f\n", x_star, y_star_L2);
    fprintf(data_file, "%.6f %.6f\n", x_star, y_star_L3);
    
    fclose(data_file);
    
    // Создаем скрипт для gnuplot
    fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
    fprintf(gnuplot_script, "set output 'lagrange_graph.png'\n");
    fprintf(gnuplot_script, "set title 'Интерполяция Лагранжа - Вариант 44' font 'Arial,14'\n");
    fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set ylabel 'y' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set grid\n");
    fprintf(gnuplot_script, "set key top left box\n");
    fprintf(gnuplot_script, "set xrange [-1.5:2.5]\n");
    fprintf(gnuplot_script, "set yrange [-3:2]\n");
    fprintf(gnuplot_script, "plot 'lagrange_data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Все исходные точки', \\\n");
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 1 with lines lw 2 lc rgb 'red' title 'L2(x) (2-я степень)', \\\n");
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'L3(x) (3-я степень)', \\\n");
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 3 with points pt 2 ps 2 lc rgb 'dark-red' title 'Узлы L2', \\\n");
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 4 with points pt 4 ps 2 lc rgb 'dark-blue' title 'Узлы L3', \\\n");
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 5 with points pt 9 ps 3 lc rgb 'green' title 'x* = %.3f'\n", x_star);
    fclose(gnuplot_script);
    
    // Запускаем gnuplot
    system("gnuplot plot_lagrange.gnu");
    printf("График сохранен в файл: lagrange_graph.png\n");
    printf("Точка x* = %.3f отмечена зеленым маркером на графике\n", x_star);
}

// Функция для вычисления разделенных разностей
void divided_differences(double x_nodes[], double y_nodes[], double diff_table[], int indices[], int n) {
    // Копируем значения y в первую колонку таблицы
    for (int i = 0; i <= n; i++) {
        diff_table[i] = y_nodes[indices[i]];
    }
    
    // Вычисляем разделенные разности
    for (int j = 1; j <= n; j++) {
        for (int i = n; i >= j; i--) {
            diff_table[i] = (diff_table[i] - diff_table[i-1]) / 
                           (x_nodes[indices[i]] - x_nodes[indices[i-j]]);
        }
    }
}

// Метод Ньютона для интерполяции
double newton_interpolation(double x, double x_nodes[], double diff_table[], int indices[], int n) {
    double result = diff_table[0];
    double product = 1.0;
    
    for (int i = 1; i <= n; i++) {
        product *= (x - x_nodes[indices[i-1]]);
        result += diff_table[i] * product;
    }
    
    return result;
}

// Оценка погрешности методом Ньютона
double newton_error_estimation(double x, double x_nodes[], double y_nodes[], int indices[], int n) {
    // Создаем таблицу разделенных разностей для n+1 узлов
    double diff_table_n[10]; // максимальный размер
    divided_differences(x_nodes, y_nodes, diff_table_n, indices, n);
    
    // Для оценки погрешности нам нужен следующий член (n+1)-го порядка
    // Берем дополнительный узел для вычисления разности (n+1)-го порядка
    if (n + 1 < 9) { // проверяем, что есть дополнительный узел
        int extended_indices[10];
        for (int i = 0; i <= n + 1; i++) {
            extended_indices[i] = (indices[0] + i) % 9; // берем последовательные узлы
        }
        
        double diff_table_extended[10];
        divided_differences(x_nodes, y_nodes, diff_table_extended, extended_indices, n + 1);
        
        // Вычисляем произведение (x - x0)(x - x1)...(x - xn)
        double product = 1.0;
        for (int i = 0; i <= n; i++) {
            product *= (x - x_nodes[indices[i]]);
        }
        
        // Оценка погрешности через разделенную разность (n+1)-го порядка
        double error_estimate = fabs(diff_table_extended[n + 1] * product);
        return error_estimate;
    }
    
    return 0.0; // если нельзя вычислить
}


void check_with_newton(double x, double x_nodes[], double y_nodes[], int indices[], int n, const char* label) {
    printf("%s\n", label);
    
    // Интерполяция Лагранжа
    double lagrange_value = lagrange(x, x_nodes, y_nodes, n, indices);
    printf("  Интерполяция Лагранжа: %.6f\n", lagrange_value);
    
    // Интерполяция Ньютона (должна дать тот же результат)
    double diff_table[10];
    divided_differences(x_nodes, y_nodes, diff_table, indices, n);
    double newton_value = newton_interpolation(x, x_nodes, diff_table, indices, n);
    printf("  Интерполяция Ньютона:  %.6f\n", newton_value);
    
    // Разница между методами
    double difference = fabs(lagrange_value - newton_value);
    printf("  Разница методов: %.2e\n", difference);
    
    // Оценка погрешности методом Ньютона
    double error_estimate1 = newton_error_estimation(x, x_nodes, y_nodes, indices, n);
    
    printf("  Оценка погрешности (метод Ньютона): %.2e\n", error_estimate1);
    
    printf("\n");
}

int main() {    
    double x[] = {-1.31, -0.86, -0.41, 0.04, 0.49, 0.94, 1.39, 1.84, 2.29};
    double y[] = {0.7814, -2.0732, -2.3769, -1.9327, 0.1364, 0.6428, -0.0213, -2.2016, -1.9837};
    
    double x_star = 1.708;
    
    printf("=== ВАРИАНТ 44 ===\n");
    printf("Точка интерполяции: x* = %.3f\n\n", x_star);
    
    printf("1.1 Многочлен Лагранжа 2-й степени (L2(x)):\n");
    int indices_L2[] = {5, 6, 7}; // точки 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 3; i++) {
        printf("(%.2f, %.4f) ", x[indices_L2[i]], y[indices_L2[i]]);
    }
    printf("\n");
    
    double L2 = lagrange(x_star, x, y, 2, indices_L2);
    printf("   L2(%.3f) = %.6f\n", x_star, L2);

    //check_with_newton(x_star, x, y, indices_L2, 2, "  Проверка методом Ньютона:");
    
    printf("\n");

    //проверка в узловой точке
    check_at_node(x, y, 6, indices_L2, 2);

    printf("1.2 Многочлен Лагранжа 2-й степени (L2(x)):\n");
    int indices_L21[] = {6, 7, 8}; // точки 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 3; i++) {
        printf("(%.2f, %.4f) ", x[indices_L21[i]], y[indices_L21[i]]);
    }
    printf("\n");
    
    double L21 = lagrange(x_star, x, y, 2, indices_L21);
    printf("   L2(%.3f) = %.6f\n", x_star, L21);

    //check_with_newton(x_star, x, y, indices_L21, 2, "  Проверка методом Ньютона:");

    printf("\n");

    //проверка в узловой точке
    check_at_node(x, y, 6, indices_L21, 2);
    
    printf("2.1 Многочлен Лагранжа 3-й степени (L3(x)):\n");
    int indices_L3[] = {4, 5, 6, 7}; // точки 0.49, 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 4; i++) {
        printf("(%.2f, %.4f) ", x[indices_L3[i]], y[indices_L3[i]]);
    }
    printf("\n");
    
    double L3 = lagrange(x_star, x, y, 3, indices_L3);
    printf("   L3(%.3f) = %.6f\n", x_star, L3);

    //check_with_newton(x_star, x, y, indices_L3, 3, "  Проверка методом Ньютона:");
    
    printf("\n");
    
    check_at_node(x, y, 6, indices_L3, 3);

    printf("2.2 Многочлен Лагранжа 3-й степени (L3(x)):\n");
    int indices_L32[] = {5, 6, 7, 8}; // точки 0.49, 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 4; i++) {
        printf("(%.2f, %.4f) ", x[indices_L32[i]], y[indices_L32[i]]);
    }
    printf("\n");
    
    double L32 = lagrange(x_star, x, y, 3, indices_L32);
    printf("   L3(%.3f) = %.6f\n", x_star, L32);

    //check_with_newton(x_star, x, y, indices_L32, 3, "  Проверка методом Ньютона:");
        
    printf("\n");
    
    check_at_node(x, y, 6, indices_L32, 3);

    printf("2.3 Многочлен Лагранжа 3-й степени (L3(x)):\n");
    int indices_L33[] = {3, 4, 5, 6}; // точки 0.49, 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 4; i++) {
        printf("(%.2f, %.4f) ", x[indices_L33[i]], y[indices_L33[i]]);
    }
    printf("\n");
    
    double L33 = lagrange(x_star, x, y, 3, indices_L33);
    printf("   L3(%.3f) = %.6f\n", x_star, L33);

    //check_with_newton(x_star, x, y, indices_L32, 3, "  Проверка методом Ньютона:");
        
    printf("\n");
    
    check_at_node(x, y, 6, indices_L33, 3);
    
    printf("3. Сравнение результатов:\n");
    printf("   L2_1(%.3f) = %.6f\n", x_star, L2);
    printf("   L2_2(%.3f) = %.6f\n", x_star, L21);
    printf("   L3_1(%.3f) = %.6f\n", x_star, L3);
    printf("   L3_2(%.3f) = %.6f\n", x_star, L32);
    printf("   L3_3(%.3f) = %.6f\n", x_star, L33);
    printf("   Разница L2 методов: %.6f\n", fabs(L2 - L21));
    printf("   Разница L3 методов: %.6f\n", fabs(L3 - L32));
    
    printf("\n4. Построение графиков...\n");
    plot_lagrange(x, y, indices_L21, indices_L33, 9, x_star);
    return 0;
}

/*
Сравнение результатов:
   L2_1(1.708) = -1.404900
   L2_2(1.708) = -1.810605
   L3_1(1.708) = -1.384516
   L3_2(1.708) = -1.635701
   L3_3(1.708) = -0.983054
*/