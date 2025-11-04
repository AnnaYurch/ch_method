#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_POINTS 9

typedef struct {
    int n;          // количество точек
    double x[MAX_POINTS];
    double y[MAX_POINTS];
    double diff[MAX_POINTS][MAX_POINTS]; // таблица разделенных разностей
} DataTable;

// вычисление разделенных разностей для произвольных узлов
void calculate_divided_differences_for_nodes(DataTable *data, int indices[], int node_count) {
    int i, j;
    
    double temp_x[MAX_POINTS], temp_y[MAX_POINTS];
    
    for (i = 0; i < node_count; i++) {
        temp_x[i] = data->x[indices[i]];
        temp_y[i] = data->y[indices[i]];
    }
    
    // инициализация первого столбца значениями функции
    for (i = 0; i < node_count; i++) {
        data->diff[i][0] = temp_y[i];
    }
    
    // вычисление разделенных разностей
    for (j = 1; j < node_count; j++) {
        for (i = 0; i < node_count - j; i++) {
            data->diff[i][j] = (data->diff[i+1][j-1] - data->diff[i][j-1]) / 
                              (temp_x[i+j] - temp_x[i]);
        }
    }
}

double newton_forward(DataTable *data, double x, int indices[], int degree) {
    if (degree >= data->n) {
        printf("Ошибка: степень слишком высокая!\n");
        return 0.0;
    }
    
    double result = data->diff[0][0]; // f(x0)
    double product = 1.0;
    int i;
    
    for (i = 1; i <= degree; i++) {
        product *= (x - data->x[indices[i-1]]); // (x-x0)(x-x1)...(x-x_{i-1})
        result += product * data->diff[0][i]; // f(x0,x1,...,xi)
    }
    
    return result;
}

double newton_backward(DataTable *data, double x, int indices[], int degree) {
    if (degree >= data->n) {
        printf("Ошибка: степень слишком высокая!\n");
        return 0.0;
    }
    
    double result = data->diff[degree][0]; // f(x_n)
    double product = 1.0;
    int i;
    
    for (i = 1; i <= degree; i++) {
        product *= (x - data->x[indices[degree - i + 1]]); // (x-x_n)(x-x_{n-1})...
        result += product * data->diff[degree - i][i]; 
    }
    
    return result;
}

// оценка погрешности интерполяции
double estimate_error(DataTable *data, double x, int indices[], int degree) {
    if (degree + 1 >= data->n) {
        return 0.0;  
    }
    
    // (x-x₀)(x-x₁)...(x-x_n)
    double product = 1.0;
    for (int i = 0; i <= degree; i++) {
        product *= (x - data->x[indices[i]]);
    }
    
    //если есть узел после последнего используемого - берем его
    //иначе берем узел перед первым используемым
    int next_index;
    if (indices[degree] + 1 < data->n) {
        next_index = indices[degree] + 1;  
    } else {
        next_index = indices[0] - 1;       
    }
    
    if (next_index < 0 || next_index >= data->n) {
        return 0.0; 
    }
    
    //создаем расширенный набор узлов
    int extended_indices[MAX_POINTS];
    for (int i = 0; i <= degree; i++) {
        extended_indices[i] = indices[i];
    }
    extended_indices[degree + 1] = next_index;
    
    //вычисляем разделенную разность (n+1)-го порядка
    calculate_divided_differences_for_nodes(data, extended_indices, degree + 2);
    
    double next_diff = data->diff[0][degree + 1];
    
    return fabs(product * next_diff);
}

void print_selected_nodes(DataTable *data, int indices[], int degree) {
    printf(" Узлы: ");
    for (int i = 0; i <= degree; i++) {
        printf("x%d=%.4f(y=%.4f) ", indices[i], 
               data->x[indices[i]], data->y[indices[i]]);
    }
    printf("\n");
}

// проверка в узловой точке
void test_at_node(DataTable *data, int node_index, int indices[], int degree, const char* method) {
    double x_node = data->x[node_index];
    double y_actual = data->y[node_index];
    double y_calculated;
    
    if (strcmp(method, "forward") == 0) {
        y_calculated = newton_forward(data, x_node, indices, degree);
    } else {
        y_calculated = newton_backward(data, x_node, indices, degree);
    }
    
    double error = fabs(y_actual - y_calculated);
    
    printf("Проверка в узле x%d=%.4f:\n", node_index, x_node);
    printf("  Фактическое значение:  y = %10.6f\n", y_actual);
    printf("  Интерполированное: P%d = %10.6f\n", degree, y_calculated);
    printf("  Погрешность:  %10.2f\n", error);
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
    printf("=== ИНТЕРПОЛЯЦИЯ НЬЮТОНА ===\n");
    printf("Точка интерполяции: x* = %.3f\n\n", x_star);
        
    printf("1.1 МНОГОЧЛЕН 2-Й СТЕПЕНИ\n");
    int nodes_P2_1[] = {5, 6, 7}; // точки 0.94, 1.39, 1.84 (ближайшие к x*)
    calculate_divided_differences_for_nodes(&data, nodes_P2_1, 3);
    
    print_selected_nodes(&data, nodes_P2_1, 2);
    double p2_1 = newton_forward(&data, x_star, nodes_P2_1, 2);
    double error_p2_1 = estimate_error(&data, x_star, nodes_P2_1, 2);
    
    printf(" P2(%.3f) = %12.8f\n", x_star, p2_1);
    printf(" Оценка погрешности: %12.8f\n\n", error_p2_1);
    
    test_at_node(&data, 6, nodes_P2_1, 2, "forward");
    
    printf("\n1.2 МНОГОЧЛЕН 2-Й СТЕПЕНИ\n");
    int nodes_P2_2[] = {6, 7, 8}; // точки 1.39, 1.84, 2.29 (конечные точки)
    calculate_divided_differences_for_nodes(&data, nodes_P2_2, 3);
    
    print_selected_nodes(&data, nodes_P2_2, 2);
    double p2_2 = newton_backward(&data, x_star, nodes_P2_2, 2);
    double error_p2_2 = estimate_error(&data, x_star, nodes_P2_2, 2);
    
    printf(" P2(%.3f) = %12.8f\n", x_star, p2_2);
    printf(" Оценка погрешности: %12.8f\n\n", error_p2_2);
    
    test_at_node(&data, 7, nodes_P2_2, 2, "backward");
    
    printf("\n 2.1 МНОГОЧЛЕН 3-Й СТЕПЕНИ\n");
    int nodes_P3_1[] = {4, 5, 6, 7}; // точки 0.49, 0.94, 1.39, 1.84
    calculate_divided_differences_for_nodes(&data, nodes_P3_1, 4);
    
    print_selected_nodes(&data, nodes_P3_1, 3);
    double p3_1 = newton_forward(&data, x_star, nodes_P3_1, 3);
    double error_p3_1 = estimate_error(&data, x_star, nodes_P3_1, 3);
    
    printf(" P3(%.3f) = %12.8f\n", x_star, p3_1);
    printf(" Оценка погрешности: %12.8f\n\n", error_p3_1);
    
    test_at_node(&data, 6, nodes_P3_1, 3, "forward");
    
    printf("\nМНОГОЧЛЕН 3-Й СТЕПЕНИ\n");
    int nodes_P3_2[] = {5, 6, 7, 8}; // точки 0.94, 1.39, 1.84, 2.29
    calculate_divided_differences_for_nodes(&data, nodes_P3_2, 4);
    
    print_selected_nodes(&data, nodes_P3_2, 3);
    double p3_2 = newton_backward(&data, x_star, nodes_P3_2, 3);
    double error_p3_2 = estimate_error(&data, x_star, nodes_P3_2, 3);
    
    printf(" P3(%.3f) = %12.8f\n", x_star, p3_2);
    printf(" Оценка погрешности: %12.8f\n\n", error_p3_2);
    
    test_at_node(&data, 8, nodes_P3_2, 3, "backward");

    printf("\n2.3 МНОГОЧЛЕН 3-Й СТЕПЕНИ\n");
    int nodes_P3_3[] = {3, 4, 5, 6}; // точки 0.04, 0.49, 0.94, 1.39
    calculate_divided_differences_for_nodes(&data, nodes_P3_3, 4);

    print_selected_nodes(&data, nodes_P3_3, 3);
    double p3_3 = newton_forward(&data, x_star, nodes_P3_3, 3);  
    double error_p3_3 = estimate_error(&data, x_star, nodes_P3_3, 3);

    printf(" P3(%.3f) = %12.8f\n", x_star, p3_3);
    printf(" Оценка погрешности: %12.8f\n\n", error_p3_3);

    test_at_node(&data, 5, nodes_P3_3, 3, "forward");  
    
    printf("\n=== СРАВНЕНИЕ РЕЗУЛЬТАТОВ ===\n");
    printf("P2 (вариант 1): %12.8f\n", p2_1);
    printf("P2 (вариант 2): %12.8f\n", p2_2);
    printf("P3 (вариант 3): %12.8f\n", p3_1);
    printf("P3 (вариант 4): %12.8f\n", p3_2);
    printf("P3 (вариант 5): %12.8f\n", p3_3);
    printf("Разница P2 методов: %12.8f\n", fabs(p2_1 - p2_2));
    printf("Разница P3 методов: %12.8f\n\n", fabs(p3_1 - p3_2));

    return 0;
}