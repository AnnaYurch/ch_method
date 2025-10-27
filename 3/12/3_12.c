#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 10  //количество точек

void gauss_solve(int n, double A[n][n+1], double x[n]) {
    //прямой ход
    for (int i = 0; i < n; i++) {
        //поиск главного элемента
        int max_row = i;
        for (int k = i+1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[max_row][i])) {
                max_row = k;
            }
        }
        
        //перестановка строк
        for (int k = i; k < n+1; k++) {
            double temp = A[i][k];
            A[i][k] = A[max_row][k];
            A[max_row][k] = temp;
        }
        
        //обнуление элементов ниже главной диагонали
        for (int k = i+1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n+1; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }
    
    //обратный ход
    for (int i = n-1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

//вычисление значения многочлена
double polynomial_value(double x, double coefs[], int degree) {
    double result = 0;
    for (int i = 0; i <= degree; i++) {
        result += coefs[i] * pow(x, i);
    }
    return result;
}

//вычисление суммы квадратов ошибок
double calculate_error(double x[], double y[], double coefs[], int degree) {
    double error = 0;
    for (int i = 0; i < N; i++) {
        double diff = polynomial_value(x[i], coefs, degree) - y[i];
        error += diff * diff;
    }
    return error;
}

void plot_graphs(double x[], double y[], double coefs1[], double coefs2[], double coefs3[], int n) {
    FILE *data_file = fopen("data.txt", "w");
    FILE *gnuplot_script = fopen("plot.gnu", "w");
    
    if (!data_file || !gnuplot_script) {
        printf("Ошибка создания файлов для графиков\n");
        return;
    }
    
    fprintf(data_file, "# ИСХОДНЫЕ ТОЧКИ\n");
    for (int i = 0; i < n; i++) {
        fprintf(data_file, "%.6f %.6f\n", x[i], y[i]);
    }
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 1-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs1, 1));
    }
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 2-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs2, 2));
    }
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 3-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs3, 3));
    }
    
    fclose(data_file);
    
    fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
    fprintf(gnuplot_script, "set output 'mnk_graph.png'\n");
    fprintf(gnuplot_script, "set title 'Метод Наименьших Квадратов - Вариант 44' font 'Arial,14'\n");
    fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set ylabel 'y' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set grid\n");
    fprintf(gnuplot_script, "set key top left box\n");
    fprintf(gnuplot_script, "plot 'data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Исходные данные', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 1 with lines lw 2 lc rgb 'red' title 'F1(x)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'F2(x)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 3 with lines lw 2 lc rgb 'green' title 'F3(x)'\n");
    fclose(gnuplot_script);
    
    system("gnuplot plot.gnu");
    printf("Исправленный график сохранен в файл: mnk_graph.png\n");
}

int main() {
    double x[N] = {-0.92, -0.13, 0.66, 1.45, 2.24, 3.03, 3.82, 4.61, 5.40, 6.19};
    double y[N] = {-1.1748, -1.5813, -1.6726, -1.2761, 0.0742, 1.1094, 1.5785, 1.4861, 0.9617, 0.1384};
    double x_star = 3.426;
    
    printf("МЕТОД НАИМЕНЬШИХ КВАДРАТОВ\n");
    printf("x* = %.3f\n\n", x_star);
    
    //вычисление сумм
    double sum_x[7] = {0};
    double sum_xy[4] = {0};
    double sum_y = 0;
    
    for (int i = 0; i < N; i++) {
        double x_power = 1.0;
        for (int j = 0; j < 7; j++) {
            sum_x[j] += x_power;
            x_power *= x[i];
        }
        
        double xy_power = 1.0;
        for (int j = 0; j < 4; j++) {
            sum_xy[j] += y[i] * xy_power;
            xy_power *= x[i];
        }
        
        sum_y += y[i];
    }
    
    printf("1. МНОГОЧЛЕН 1-й СТЕПЕНИ: F1(x) = a0 + a1*x\n");
    double A1[2][3] = {
        {sum_x[0], sum_x[1], sum_y},      // 10·a0 + 26.35·a1 = -0.3565
        {sum_x[1], sum_x[2], sum_xy[1]}   // 26.35·a0 + 120.92·a1 = 20.7905
    };
    
    double coefs1[2];
    gauss_solve(2, A1, coefs1);
    
    printf("Коэффициенты: a0 = %.6f, a1 = %.6f\n", coefs1[0], coefs1[1]);
    printf("F1(x) = %.6f + %.6f*x\n", coefs1[0], coefs1[1]);
    
    double error1 = calculate_error(x, y, coefs1, 1);
    printf("Сумма квадратов ошибок: Φ1 = %.6f\n", error1);
    
    double F1_xstar = polynomial_value(x_star, coefs1, 1);
    printf("F1(%.3f) = %.6f\n\n", x_star, F1_xstar);
    
    printf("2. МНОГОЧЛЕН 2-й СТЕПЕНИ: F2(x) = a0 + a1*x + a2*x^2\n");
    double A2[3][4] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_y},
        {sum_x[1], sum_x[2], sum_x[3], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_xy[2]}
    };
    
    double coefs2[3];
    gauss_solve(3, A2, coefs2);
    
    printf("Коэффициенты: a0 = %.6f, a1 = %.6f, a2 = %.6f\n", coefs2[0], coefs2[1], coefs2[2]);
    printf("F2(x) = %.6f + %.6f*x + %.6f*x^2\n", coefs2[0], coefs2[1], coefs2[2]);
    
    double error2 = calculate_error(x, y, coefs2, 2);
    printf("Сумма квадратов ошибок: Φ2 = %.6f\n", error2);
    
    double F2_xstar = polynomial_value(x_star, coefs2, 2);
    printf("F2(%.3f) = %.6f\n\n", x_star, F2_xstar);
    
    printf("3. МНОГОЧЛЕН 3-й СТЕПЕНИ: F3(x) = a0 + a1*x + a2*x^2 + a3*x^3\n");
    double A3[4][5] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_x[3], sum_y},
        {sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_xy[2]},
        {sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_xy[3]}
    };
    
    double coefs3[4];
    gauss_solve(4, A3, coefs3);
    
    printf("Коэффициенты: a0 = %.6f, a1 = %.6f, a2 = %.6f, a3 = %.6f\n", coefs3[0], coefs3[1], coefs3[2], coefs3[3]);
    printf("F3(x) = %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3\n", coefs3[0], coefs3[1], coefs3[2], coefs3[3]);
    
    double error3 = calculate_error(x, y, coefs3, 3);
    printf("Сумма квадратов ошибок: Φ3 = %.6f\n", error3);
    
    double F3_xstar = polynomial_value(x_star, coefs3, 3);
    printf("F3(%.3f) = %.6f\n\n", x_star, F3_xstar);
    
    printf("СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ:\n");
    printf("Степень | Сумма квадратов ошибок | F(x*)\n");
    printf("--------|------------------------|---------\n");
    printf("   1    |       %12.6f     | %8.4f\n", error1, F1_xstar);
    printf("   2    |       %12.6f     | %8.4f\n", error2, F2_xstar);
    printf("   3    |       %12.6f     | %8.4f\n", error3, F3_xstar);
    
    printf("\nСтроим графики...\n");
    plot_graphs(x, y, coefs1, coefs2, coefs3, N);
    
    return 0;
}