#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define N 5
#define EPSILON 0.0001
#define MAX_ITER 1000

typedef struct {
    double real;
    double imag;
} Complex;

void print_matrix(double A[N][N], const char* name) {
    printf("%s:\n", name);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%10.6f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_complex(Complex c, const char* name) {
    if (fabs(c.imag) < 1e-10) {
        printf("%s = %.6f\n", name, c.real);
    } else if (c.imag > 0) {
        printf("%s = %.6f + %.6fi\n", name, c.real, c.imag);
    } else {
        printf("%s = %.6f - %.6fi\n", name, c.real, -c.imag);
    }
}

void copy_matrix(double src[N][N], double dest[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

void matrix_mult(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void transpose(double A[N][N], double AT[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            AT[j][i] = A[i][j];
        }
    }
}

double vector_norm(double v[N]) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

//алгоритм Хаусхолдера
void householder_qr(double A[N][N], double Q[N][N], double R[N][N]) {
    double A_temp[N][N];
    copy_matrix(A, A_temp);
    
    //инициализация Q как единичной матрицы
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    for (int k = 0; k < N - 1; k++) {
        double x[N], v[N], v_norm;
        
        //формирование вектора x
        for (int i = 0; i < N; i++) {
            if (i < k) { //элементы выше текущего столбца не трогаем
                x[i] = 0;
            } else {
                x[i] = A_temp[i][k];
            }
        }
        
        double x_norm = vector_norm(x);
        
        //формирование вектора v
        for (int i = 0; i < N; i++) {
            if (i < k) {
                v[i] = 0;
            } else if (i == k) {
                v[i] = x[i] + copysign(x_norm, x[i]); //copysign(5.0, -3.0) // = -5.0 
            } else {
                v[i] = x[i];
            }
        }
        
        v_norm = vector_norm(v);
        
        //нормализация v
        if (v_norm > 1e-12) {
            for (int i = k; i < N; i++) {
                v[i] /= v_norm;
            }
        }
        
        //построение H
        //единичная матрица
        double H[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    H[i][j] = 1.0;
                } else {
                    H[i][j] = 0.0;
                }
            }
        }

        // H = I - 2·v·vᵀ
        for (int i = k; i < N; i++) {
            for (int j = k; j < N; j++) {
                H[i][j] -= 2 * v[i] * v[j];
            }
        }
        
        // A_temp1 = H * A_temp0 (для получения нулей под диагональю)
        double temp[N][N];
        matrix_mult(H, A_temp, temp);
        copy_matrix(temp, A_temp);
        
        //обновление Q: Q = Q * H (H1*H2*H3*....)
        matrix_mult(Q, H, temp);
        copy_matrix(temp, Q);
    }
    
    copy_matrix(A_temp, R);
}

//нахождения собственных значений
void qr_algorithm(double A[N][N], Complex eigenvalues[N]) {
    double A_k[N][N], Q[N][N], R[N][N], A_next[N][N];
    double QT[N][N];
    
    copy_matrix(A, A_k);
    
    for (int iter = 0; iter < MAX_ITER; iter++) {
        householder_qr(A_k, Q, R);
        
        // A_{k+1} = R * Q
        matrix_mult(R, Q, A_next);
        
        //проверка сходимости (по поддиагональным элементам)
        bool converged = true;
        for (int m = 0; m < N - 1; m++) {
            double sum_squares = 0.0;
            //вычисляем сумму(a_im)²
            for (int i = m + 1; i < N; i++) {
                sum_squares += A_next[i][m] * A_next[i][m];
            }
            double norm = sqrt(sum_squares);
            
            if (norm >= EPSILON) {
                converged = false;
                break;
            }
        } //критерий сходимости - матрица становится верхней треугольной
        
        copy_matrix(A_next, A_k);
        
        if (converged) {
            printf("Сходимость достигнута на итерации %d\n\n", iter + 1);
            break;
        }
    }

    printf("Последняя квазитреугольная матрица A_final:\n");
    printf("==========================================\n");
    print_matrix(A_k, "A_final");

    //извлечение собственных значений из квазитреугольной матрицы
    int i = 0;
    while (i < N) {
        if (i == N - 1 || fabs(A_k[i+1][i]) < EPSILON) {
            //вещественное собственное значение
            eigenvalues[i].real = A_k[i][i];
            eigenvalues[i].imag = 0;
            i++;
        } else {
            //комплексно-сопряженная пара
            double a = A_k[i][i];
            double b = A_k[i][i+1];
            double c = A_k[i+1][i];
            double d = A_k[i+1][i+1];
            
            // λ² - (a+d)λ + (ad - bc) = 0
            double trace = a + d;
            double det = a * d - b * c;
            double discriminant = trace * trace - 4 * det;
            
            if (discriminant >= 0) {
                //вещественные корни
                eigenvalues[i].real = (trace + sqrt(discriminant)) / 2;
                eigenvalues[i].imag = 0;
                eigenvalues[i+1].real = (trace - sqrt(discriminant)) / 2;
                eigenvalues[i+1].imag = 0;
            } else {
                //комплексно-сопряженные корни
                //λ₁ = α + βi
                //λ₂ = α - βi
                eigenvalues[i].real = trace / 2;
                eigenvalues[i].imag = sqrt(-discriminant) / 2;
                eigenvalues[i+1].real = trace / 2;
                eigenvalues[i+1].imag = -sqrt(-discriminant) / 2;
            }
            i += 2;
        }
    }
}

int main() {
    double A[N][N] = {
        {4, -5, -4, 7, 5},
        {1, 1, 0, 1, 1},
        {-2, 3, 4, 8, 5},
        {2, 5, -4, 0, 3},
        {6, 9, -5, 1, -1}
    };
    
    printf("QR-алгоритм для нахождения собственных значений\n");
    printf("Точность: %f\n", EPSILON);
    printf("===============================================\n\n");
    
    print_matrix(A, "Исходная матрица A");
    
    //проверка QR-разложения на первой итерации
    double Q[N][N], R[N][N];
    householder_qr(A, Q, R);
    
    print_matrix(Q, "Матрица Q после QR-разложения");
    print_matrix(R, "Матрица R после QR-разложения");
    
    //проверка: A = Q * R
    double QR[N][N];
    matrix_mult(Q, R, QR);
    print_matrix(QR, "Проверка: Q * R");
    
    //нахождение собственных значений
    Complex eigenvalues[N];
    qr_algorithm(A, eigenvalues);
    
    printf("Собственные значения:\n");
    printf("=====================\n");
    for (int i = 0; i < N; i++) {
        char name[20];
        sprintf(name, "λ%d", i + 1);
        print_complex(eigenvalues[i], name);
    }
    
    return 0;
}