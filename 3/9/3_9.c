#include <stdio.h>
#include <math.h>

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
    
    return 0;
}