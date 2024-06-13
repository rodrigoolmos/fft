#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define DIM 1024*64

typedef struct {
    float real;
    float imag;
} Complex;

void FFT(Complex *x, int N) {
    int i, j, m, len;
    float angle;
    Complex temp, wlen, w, u, v, next_w;

    j = 0;
    for (i = 0; i < N; ++i) { // bit reversal
        if (i < j) {
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
        m = N >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    for (len = 2; len <= N; len <<= 1) {
        angle = -2 * PI / len;
        wlen.real = cos(angle);
        wlen.imag = sin(angle);
        for (i = 0; i < N; i += len) {
            w.real = 1.0;
            w.imag = 0.0;
            for (j = 0; j < len / 2; ++j) {
                u = x[i + j];
                v.real = x[i + j + len / 2].real * w.real - x[i + j + len / 2].imag * w.imag;
                v.imag = x[i + j + len / 2].real * w.imag + x[i + j + len / 2].imag * w.real;
                x[i + j].real = u.real + v.real;
                x[i + j].imag = u.imag + v.imag;
                x[i + j + len / 2].real = u.real - v.real;
                x[i + j + len / 2].imag = u.imag - v.imag;
                next_w.real = w.real * wlen.real - w.imag * wlen.imag;
                next_w.imag = w.real * wlen.imag + w.imag * wlen.real;
                w = next_w;
            }
        }
    }
}


void DFT(Complex *input, Complex *output, int N) {
    for (int k = 0; k < N; k++) {  // Para cada punto de salida
        output[k].real = 0;
        output[k].imag = 0;
        for (int n = 0; n < N; n++) {  // Para cada punto de entrada
            float angle = 2 * PI * k * n / N;  // Cambio de signo en el ángulo
            output[k].real += input[n].real * cos(angle) - input[n].imag * sin(angle);
            output[k].imag += input[n].imag * cos(angle) + input[n].real * sin(angle);
        }
    }
}

void printComplexArray(Complex *array, int N) {
    for (int i = 0; i < N; i++) {
        printf("(%f, %f)\n", array[i].real, array[i].imag);
    }
}

int main() {
    int N = DIM;  // Número de puntos

    // Array de entrada de ejemplo (4 puntos)
    Complex input[DIM];

    for (int i = 0; i < N; i++) {
        input[i].real = (float)rand();
        input[i].imag = (float)rand();
    }

    // Array de salida
    Complex output[N];

    // Calcular la DFT
    clock_t start = clock();
    DFT(input, output, N);
    clock_t end = clock();

    printf("Tiempo de ejecución: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    FFT(input, N);
    end = clock();

    printf("Tiempo de ejecución: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);


    if(memcmp(input, output, sizeof(Complex)* N)){
        printf("OK\n");
    }else{
        printf("ERR\n");
    }

    // Imprimir la DFT

    return 0;
}
