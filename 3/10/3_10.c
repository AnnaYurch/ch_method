#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 9

typedef struct {
    int n;          //количество точек
    double x[MAX_POINTS];
    double y[MAX_POINTS];
    double diff[MAX_POINTS][MAX_POINTS]; //таблица разделенных разностей
} DataTable;

//вычисления разделенных разностей
void calculate_divided_differences(DataTable *data) {
    int i, j;
    
    //инициализация первого столбца значениями функции
    for (i = 0; i < data->n; i++) {
        data->diff[i][0] = data->y[i];
    }
    
    //вычисление разделенных разностей
    for (j = 1; j < data->n; j++) {
        for (i = 0; i < data->n - j; i++) {
            data->diff[i][j] = (data->diff[i+1][j-1] - data->diff[i][j-1]) / 
                              (data->x[i+j] - data->x[i]);
        }
    }
}

// Многочлен Ньютона с явным указанием начального узла и степени
double newton_custom(DataTable *data, double x, int start_index, int degree) {
    if (start_index < 0 || start_index + degree >= data->n) {
        printf("Ошибка: неверные индексы узлов!\n");
        return 0.0;
    }
    
    double result = data->diff[start_index][0];
    double product = 1.0;
    int i;
    
    for (i = 1; i <= degree; i++) {
        product *= (x - data->x[start_index + i - 1]);
        result += product * data->diff[start_index][i];
    }
    
    return result;
}

// Оценка погрешности для кастомных узлов
double estimate_custom_error(DataTable *data, double x, int start_index, int degree) {
    if (start_index < 0 || start_index + degree >= data->n) {
        return 0.0;
    }
    
    double product = 1.0;
    int i;
    
    for (i = 0; i <= degree; i++) {
        product *= (x - data->x[start_index + i]);
    }
    
    // Используем следующую разделенную разность для оценки погрешности
    double next_diff;
    if (start_index + degree + 1 < data->n) {
        next_diff = data->diff[start_index][degree + 1];
    } else {
        next_diff = data->diff[data->n - degree - 2][degree + 1];
    }
    
    return fabs(product * next_diff);
}

// Вывод информации об используемых узлах
void print_used_nodes(DataTable *data, int start_index, int degree) {
    printf("Используемые узлы: ");
    for (int i = 0; i <= degree; i++) {
        printf("x%d=%.4f(y=%.4f) ", start_index + i, 
               data->x[start_index + i], data->y[start_index + i]);
    }
    printf("\n");
}

// ПРОВЕРКА в узловой точке (НОВАЯ ФУНКЦИЯ)
void test_at_node(DataTable *data, int node_index, int start_index, int degree) {
    double x_node = data->x[node_index];
    double y_actual = data->y[node_index];
    double y_calculated = newton_custom(data, x_node, start_index, degree);
    double error = fabs(y_actual - y_calculated);
    
    printf("Проверка в узле x%d=%.4f:\n", node_index, x_node);
    printf(" Фактическое:  y = %10.6f\n", y_actual);
    printf(" Вычисленное: P%d = %10.6f\n", degree, y_calculated);
    printf(" Погрешность:     %10.2e\n", error);
    printf("\n");

    if (error < 1e-10) {
        printf("Многочлен точно проходит через узел\n");
    } else {
        printf("Замечание: погрешность значительная\n");
    }
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
    
    //вычисление разделенных разностей
    calculate_divided_differences(&data);
    
    printf("Точка интерполяции: x* = %.3f\n\n", x_star);
    
    printf("=== УЗЛЫ: 5,6,7 ===\n");
    
    int manual_start_1 = 5; 
    double manual_p2_1 = newton_custom(&data, x_star, manual_start_1, 2);
    double error_p2_1 = estimate_custom_error(&data, x_star, manual_start_1, 2);
    
    print_used_nodes(&data, manual_start_1, 2);
    printf(" P2(%.3f) = %12.8f\n", x_star, manual_p2_1);
    printf(" Оценка погрешности: %10.2e\n", error_p2_1);
    
    printf("=== УЗЛЫ: 4,5,6,7 ===\n");
    
    int manual_start_2 = 4;
    double manual_p3_1 = newton_custom(&data, x_star, manual_start_2, 3);
    double error_p3_1 = estimate_custom_error(&data, x_star, manual_start_2, 3);
    
    print_used_nodes(&data, manual_start_2, 3);
    printf(" P3(%.3f) = %12.8f\n", x_star, manual_p3_1);
    printf(" Оценка погрешности: %10.2e\n", error_p3_1);
    
    printf("\n");

    test_at_node(&data, 5, manual_start_2, 3);

    return 0;
}