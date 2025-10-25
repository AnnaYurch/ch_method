#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 9

typedef struct {
    int n;          // количество точек
    double x[MAX_POINTS];
    double y[MAX_POINTS];
    double diff[MAX_POINTS][MAX_POINTS]; // таблица разделенных разностей
} DataTable;

// вычисления разделенных разностей
void calculate_divided_differences(DataTable *data) {
    int i, j;
    
    // инициализация первого столбца значениями функции
    for (i = 0; i < data->n; i++) {
        data->diff[i][0] = data->y[i];
    }
    
    // вычисление разделенных разностей
    for (j = 1; j < data->n; j++) {
        for (i = 0; i < data->n - j; i++) {
            data->diff[i][j] = (data->diff[i+1][j-1] - data->diff[i][j-1]) / 
                              (data->x[i+j] - data->x[i]);
        }
    }
}

double newton_backward(DataTable *data, double x, int degree) {
    if (degree >= data->n) {
        printf("Ошибка: степень слишком высокая!\n");
        return 0.0;
    }
    
    double result = data->diff[data->n-1][0]; 
    double product = 1.0;
    int i;
    
    for (i = 1; i <= degree; i++) {
        product *= (x - data->x[data->n - i]);  
        result += product * data->diff[data->n - i - 1][i]; 
    }
    return result;
}

//оценка погрешности
double estimate_backward_error(DataTable *data, double x, int degree) {
    if (degree >= data->n) {
        return 0.0;
    }
    
    double product = 1.0;
    int i;
    
    // (x-x8)(x-x7)(x-x6)
    for (i = 0; i <= degree; i++) {
        product *= (x - data->x[data->n - 1 - i]);
    }
    
    double next_diff;
    if (data->n - degree - 2 >= 0) {
        next_diff = data->diff[data->n - degree - 2][degree + 1];
    } else {
        next_diff = data->diff[0][degree + 1];
    }
    
    return fabs(product * next_diff);
}

void print_backward_nodes(DataTable *data, int degree) {
    printf("Используемые узлы: ");
    for (int i = 0; i <= degree; i++) {
        int node_index = data->n - 1 - i;
        printf("x%d=%.4f(y=%.4f) ", node_index, 
               data->x[node_index], data->y[node_index]);
    }
    printf("\n");
}

// ПРОВЕРКА в узловой точке
void test_backward_at_node(DataTable *data, int node_index, int degree) {
    double x_node = data->x[node_index];
    double y_actual = data->y[node_index];
    double y_calculated = newton_backward(data, x_node, degree);
    double error = fabs(y_actual - y_calculated);
    
    printf("\n");

    printf("Проверка в узле x%d=%.4f:\n", node_index, x_node);
    printf(" Фактическое:  y = %10.6f\n", y_actual);
    printf(" Вычисленное: P%d = %10.6f\n", degree, y_calculated);
    printf(" Погрешность:     %10.2e\n", error);
    
    if (error < 1e-10) {
        printf("Многочлен точно проходит через узел\n");
    } else {
        printf("Замечание: погрешность значительная\n");
    }
    printf("\n");
}

int main() {
    DataTable data;
    double x_star = 1.708;
    
    data.n = 9;
    double x_values[] = {-1.31, -0.86, -0.41, 0.04, 0.49, 0.94, 1.39, 1.84, 2.29};
    double y_values[] = {0.7814, -2.0732, -2.3769, -1.9327, 0.1364, 0.6428, -0.0213, -2.2016, -1.9837};
    
    for (int i = 0; i < data.n; i++) {
        data.x[i] = x_values[i];
        data.y[i] = y_values[i];
    }
    
    printf("=== ВАРИАНТ 44 ===\n");
    printf("=== ВТОРАЯ ФОРМУЛА НЬЮТОНА ===\n");
    
    // вычисление разделенных разностей
    calculate_divided_differences(&data);
    
    printf("Точка интерполяции: x* = %.3f\n\n", x_star);
    
    // МНОГОЧЛЕН 2-Й СТЕПЕНИ 
    printf("=== МНОГОЧЛЕН 2-Й СТЕПЕНИ ===\n");
    double p2 = newton_backward(&data, x_star, 2);
    double error_p2 = estimate_backward_error(&data, x_star, 2);
    
    print_backward_nodes(&data, 2);  
    printf(" P2(%.3f) = %12.8f\n", x_star, p2);
    printf(" Оценка погрешности: %12.8f\n", error_p2);
    

    test_backward_at_node(&data, 7, 2);
    
    // МНОГОЧЛЕН 3-Й СТЕПЕНИ
    printf("=== МНОГОЧЛЕН 3-Й СТЕПЕНИ ===\n");
    double p3 = newton_backward(&data, x_star, 3);
    double error_p3 = estimate_backward_error(&data, x_star, 3);
    
    print_backward_nodes(&data, 3);
    printf(" P3(%.3f) = %12.8f\n", x_star, p3);
    printf(" Оценка погрешности: %12.8f\n", error_p3);
    
    test_backward_at_node(&data, 8, 3);  

    return 0;
}