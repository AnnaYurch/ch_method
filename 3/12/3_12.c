#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//показать н аграфики ошибку с 4 степенью
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

void plot_graphs(double x[], double y[], double coefs1[], double coefs2[], double coefs3[], double coefs4[], double coefs5[], double coefs6[], int n) {
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
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 4-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs4, 4));
    }
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 5-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs5, 5));
    }
    fprintf(data_file, "\n\n");
    
    fprintf(data_file, "# МНОГОЧЛЕН 6-й СТЕПЕНИ\n");
    for (double xi = -1.0; xi <= 6.5; xi += 0.1) {
        fprintf(data_file, "%.6f %.6f\n", xi, polynomial_value(xi, coefs6, 6));
    }
    
    fclose(data_file);
    
    fprintf(gnuplot_script, "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n");
    fprintf(gnuplot_script, "set output 'mnk_graph.png'\n");
    fprintf(gnuplot_script, "set title 'Метод Наименьших Квадратов - Аппроксимация многочленами 1-6 степени' font 'Arial,14'\n");
    fprintf(gnuplot_script, "set xlabel 'x' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set ylabel 'y' font 'Arial,12'\n");
    fprintf(gnuplot_script, "set grid\n");
    fprintf(gnuplot_script, "set key top left box\n");
    fprintf(gnuplot_script, "plot 'data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Исходные данные', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 1 with lines lw 2 lc rgb 'red' title 'F1(x) (1 степень)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'F2(x) (2 степень)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 3 with lines lw 2 lc rgb 'green' title 'F3(x) (3 степень)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 4 with lines lw 2 lc rgb 'purple' title 'F4(x) (4 степень)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 5 with lines lw 2 lc rgb 'orange' title 'F5(x) (5 степень)', \\\n");
    fprintf(gnuplot_script, "     'data.txt' index 6 with lines lw 2 lc rgb 'brown' title 'F6(x) (6 степень)'\n");
    fclose(gnuplot_script);
    
    system("gnuplot plot.gnu");
    printf("График сохранен в файл: mnk_graph.png\n");
}

int main() {
    double x[N] = {-0.92, -0.13, 0.66, 1.45, 2.24, 3.03, 3.82, 4.61, 5.40, 6.19};
    double y[N] = {-1.1748, -1.5813, -1.6726, -1.2761, 0.0742, 1.1094, 1.5785, 1.4861, 0.9617, 0.1384};
    double x_star = 3.426;
    
    printf("МЕТОД НАИМЕНЬШИХ КВАДРАТОВ\n");
    printf("x* = %.3f\n\n", x_star);
    
    //вычисление сумм
    double sum_x[13] = {0};  
    double sum_xy[7] = {0};  
    double sum_y = 0;

    for (int i = 0; i < N; i++) {
        double x_power = 1.0;
        for (int j = 0; j < 13; j++) {  
            sum_x[j] += x_power;
            x_power *= x[i];
        }
        
        double xy_power = 1.0;
        for (int j = 0; j < 7; j++) {   
            sum_xy[j] += y[i] * xy_power;
            xy_power *= x[i];
        }
        
        sum_y += y[i];
    }
    
    printf("1. МНОГОЧЛЕН 1-й СТЕПЕНИ: F1(x) = a0 + a1*x\n");
    double A1[2][3] = {
        {sum_x[0], sum_x[1], sum_y},      
        {sum_x[1], sum_x[2], sum_xy[1]}   
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

    printf("4. МНОГОЧЛЕН 4-й СТЕПЕНИ: F4(x) = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4\n");
    double A4[5][6] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_y},
        {sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_xy[2]},
        {sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_xy[3]},
        {sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_xy[4]}
    };

    double coefs4[5];
    gauss_solve(5, A4, coefs4);

    printf("Коэффициенты: a0 = %.6f, a1 = %.6f, a2 = %.6f, a3 = %.6f, a4 = %.6f\n", 
        coefs4[0], coefs4[1], coefs4[2], coefs4[3], coefs4[4]);
    printf("F4(x) = %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3 + %.6f*x^4\n", 
        coefs4[0], coefs4[1], coefs4[2], coefs4[3], coefs4[4]);

    double error4 = calculate_error(x, y, coefs4, 4);
    printf("Сумма квадратов ошибок: Φ4 = %.6f\n", error4);

    double F4_xstar = polynomial_value(x_star, coefs4, 4);
    printf("F4(%.3f) = %.6f\n\n", x_star, F4_xstar);

    printf("5. МНОГОЧЛЕН 5-й СТЕПЕНИ: F5(x) = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5\n");
    double A5[6][7] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_y},
        {sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_xy[2]},
        {sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_xy[3]},
        {sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_xy[4]},
        {sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_x[10], sum_xy[5]}
    };

    double coefs5[6];
    gauss_solve(6, A5, coefs5);

    printf("Коэффициенты: a0 = %.6f, a1 = %.6f, a2 = %.6f, a3 = %.6f, a4 = %.6f, a5 = %.6f\n", 
        coefs5[0], coefs5[1], coefs5[2], coefs5[3], coefs5[4], coefs5[5]);
    printf("F5(x) = %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3 + %.6f*x^4 + %.6f*x^5\n", 
        coefs5[0], coefs5[1], coefs5[2], coefs5[3], coefs5[4], coefs5[5]);

    double error5 = calculate_error(x, y, coefs5, 5);
    printf("Сумма квадратов ошибок: Φ5 = %.6f\n", error5);

    double F5_xstar = polynomial_value(x_star, coefs5, 5);
    printf("F5(%.3f) = %.6f\n\n", x_star, F5_xstar);

    printf("6. МНОГОЧЛЕН 6-й СТЕПЕНИ: F6(x) = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6\n");
    double A6[7][8] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_y},
        {sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_xy[2]},
        {sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_xy[3]},
        {sum_x[4], sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_x[10], sum_xy[4]},
        {sum_x[5], sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_x[10], sum_x[11], sum_xy[5]},
        {sum_x[6], sum_x[7], sum_x[8], sum_x[9], sum_x[10], sum_x[11], sum_x[12], sum_xy[6]}
    };

    double coefs6[7];
    gauss_solve(7, A6, coefs6);

    printf("Коэффициенты:\n");
    printf("a0 = %.6f, a1 = %.6f, a2 = %.6f\n", coefs6[0], coefs6[1], coefs6[2]);
    printf("a3 = %.6f, a4 = %.6f, a5 = %.6f, a6 = %.6f\n", coefs6[3], coefs6[4], coefs6[5], coefs6[6]);
    printf("F6(x) = %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3 + %.6f*x^4 + %.6f*x^5 + %.6f*x^6\n", 
        coefs6[0], coefs6[1], coefs6[2], coefs6[3], coefs6[4], coefs6[5], coefs6[6]);

    double error6 = calculate_error(x, y, coefs6, 6);
    printf("Сумма квадратов ошибок: Φ6 = %.6f\n", error6);

    double F6_xstar = polynomial_value(x_star, coefs6, 6);
    printf("F6(%.3f) = %.6f\n\n", x_star, F6_xstar);
    
    printf("СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ:\n");
    printf("Степень | Сумма квадратов ошибок | F(x*)\n");
    printf("--------|------------------------|---------\n");
    printf("   1    |       %12.6f     | %8.4f\n", error1, F1_xstar);
    printf("   2    |       %12.6f     | %8.4f\n", error2, F2_xstar);
    printf("   3    |       %12.6f     | %8.4f\n", error3, F3_xstar);
    printf("   4    |       %12.6f     | %8.4f\n", error4, F4_xstar);
    printf("   5    |       %12.6f     | %8.4f\n", error5, F5_xstar);
    printf("   6    |       %12.6f     | %8.4f\n", error6, F6_xstar);

    printf("\nСтроим графики...\n");
    plot_graphs(x, y, coefs1, coefs2, coefs3, coefs4, coefs5, coefs6, N);
    
    return 0;
}