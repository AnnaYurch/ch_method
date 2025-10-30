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

void plot_lagrange(double x_nodes[], double y_nodes[], int indices_L2[], int indices_L3[], int n_total) {
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
    fprintf(gnuplot_script, "     'lagrange_data.txt' index 4 with points pt 4 ps 2 lc rgb 'dark-blue' title 'Узлы L3'\n");
    fclose(gnuplot_script);
    
    // Запускаем gnuplot
    system("gnuplot plot_lagrange.gnu");
    printf("График сохранен в файл: lagrange_graph.png\n");
}

int main() {    
    double x[] = {-1.31, -0.86, -0.41, 0.04, 0.49, 0.94, 1.39, 1.84, 2.29};
    double y[] = {0.7814, -2.0732, -2.3769, -1.9327, 0.1364, 0.6428, -0.0213, -2.2016, -1.9837};
    
    double x_star = 1.708;
    
    printf("=== ВАРИАНТ 44 ===\n");
    printf("Точка интерполяции: x* = %.3f\n\n", x_star);
    
    printf("1. Многочлен Лагранжа 2-й степени (L2(x)):\n");
    int indices_L2[] = {5, 6, 7}; // точки 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 3; i++) {
        printf("(%.2f, %.4f) ", x[indices_L2[i]], y[indices_L2[i]]);
    }
    printf("\n");
    
    double L2 = lagrange(x_star, x, y, 2, indices_L2);
    printf("   L2(%.3f) = %.6f\n", x_star, L2);
    
    printf("\n");

    //проверка в узловой точке
    check_at_node(x, y, 6, indices_L2, 2);
    
    printf("2. Многочлен Лагранжа 3-й степени (L3(x)):\n");
    int indices_L3[] = {4, 5, 6, 7}; // точки 0.49, 0.94, 1.39, 1.84
    printf("   Узлы: ");
    for (int i = 0; i < 4; i++) {
        printf("(%.2f, %.4f) ", x[indices_L3[i]], y[indices_L3[i]]);
    }
    printf("\n");
    
    double L3 = lagrange(x_star, x, y, 3, indices_L3);
    printf("   L3(%.3f) = %.6f\n", x_star, L3);
    
    printf("\n");
    
    check_at_node(x, y, 6, indices_L3, 3);
    
    printf("3. Сравнение результатов:\n");
    printf("   L2(%.3f) = %.6f\n", x_star, L2);
    printf("   L3(%.3f) = %.6f\n", x_star, L3);
    printf("   Разница: %.6f\n", fabs(L2 - L3));
    
    printf("\n4. Построение графиков...\n");
    plot_lagrange(x, y, indices_L2, indices_L3, 9);

    return 0;
}