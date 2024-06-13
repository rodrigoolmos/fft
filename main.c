#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define DIM 1024
#define EPSILON 0.2

typedef struct {
    float real;
    float imag;
} Complex;

void FFT(Complex *x, int N) {
    int i, j, m, len;
    float angle;
    Complex temp, wlen, w, u, v, next_w;

    // Reordenamiento por bit reversal
    j = 0;
    for (i = 0; i < N; ++i) {
        if (i < j) {
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
        for (m = N >> 1; m >= 1 && j >= m; m >>= 1) {
            j -= m;
        }
        j += m;
    }

    // Algoritmo de la mariposa
    for (len = 2; len <= N; len <<= 1) {
        angle = 2 * PI / len;
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
    for (int k = 0; k < N; k++) {
        output[k].real = 0;
        output[k].imag = 0;
        for (int n = 0; n < N; n++) {
            float angle = 2 * PI * k * n / N;
            output[k].real += input[n].real * cos(angle) - input[n].imag * sin(angle);
            output[k].imag += input[n].imag * cos(angle) + input[n].real * sin(angle);
        }
    }
}

int main() {
    int N = DIM;  // Número de puntos
    int error = 0;
    double real_des, real_max;
    double imag_des, imag_max;

    // Array de entrada de ejemplo (DIM puntos)
    Complex input[DIM];
    for (int i = 0; i < N; i++) {
        input[i].real = (float)rand() * 0.456;
        input[i].imag = (float)rand() * 0.345;
    }

    Complex outputDFT[DIM];
    Complex outputFFT[DIM];

    // Copia de entrada para FFT
    memcpy(outputFFT, input, sizeof(Complex) * DIM);

    clock_t start = clock();
    DFT(input, outputDFT, N);
    clock_t end = clock();
    printf("Tiempo de ejecución de la DFT: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    FFT(outputFFT, N);  // La función FFT modifica outputFFT en su lugar
    end = clock();
    printf("Tiempo de ejecución de la FFT: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);

    real_max = 0;
    imag_max = 0;
    for (int i = 0; i < N; i++) {
        real_des = fabs((outputFFT[i].real - outputDFT[i].real) / outputDFT[i].real);
        imag_des = fabs((outputFFT[i].imag - outputDFT[i].imag) / outputDFT[i].imag);
        if(real_des > EPSILON || imag_des > EPSILON)
            error = 1;

        if(real_max < real_des)
            real_max = real_des;
        if (imag_max < imag_des)
            imag_max = imag_des;        
    }

    if(error){
        printf("ERR\n");
    } else {
        printf("OK\n");
    }

    printf("MAX DIF REAL %f MAX DIF IMAG %f\n", real_max, imag_max);

    return 0;
}
