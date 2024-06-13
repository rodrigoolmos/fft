#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI M_PI
#define DIM 512
#define EPSILON 0.002

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
        angle = - 2 * PI / len;
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
            float angle = - 2 * PI * k * n / N;
            output[k].real += input[n].real * cos(angle) - input[n].imag * sin(angle);
            output[k].imag += input[n].imag * cos(angle) + input[n].real * sin(angle);
        }
    }
}

int main() {
    int N = DIM;  // Número de puntos
    float dif, total_dif;

    // Array de entrada de ejemplo (DIM puntos)
    Complex input[DIM];
    for (int i = 0; i < N; i++) {
        input[i].real = 5*sin(2*PI*i/12) + 3*sin(2*PI*i/43) + 2*sin(2*PI*i/54);
        input[i].imag = 3*sin(2*PI*i/61) + 6*sin(2*PI*i/23) + 8*sin(2*PI*i/14);
    }

    Complex outputDFT[DIM];
    Complex outputFFT[DIM];

    // Copia de entrada para FFT
    memcpy(outputFFT, input, sizeof(Complex) * DIM);

    clock_t start = clock();
    DFT(input, outputDFT, N);
    clock_t end = clock();
    //printf("Tiempo de ejecución de la DFT: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    FFT(outputFFT, N);  // La función FFT modifica outputFFT en su lugar
    end = clock();
    //printf("Tiempo de ejecución de la FFT: %f segundos\n", (double) (end - start) / CLOCKS_PER_SEC);

    total_dif = 0;
    for (int i = 0; i < N; i++) {

        dif = fabs(sqrtf((outputFFT[i].real * outputFFT[i].real) + (outputFFT[i].imag * outputFFT[i].imag)) -
              sqrtf((outputDFT[i].real * outputDFT[i].real) + (outputDFT[i].imag * outputDFT[i].imag)));

        total_dif += dif;
    }

    printf("Input = single([");
    for (int e = 0; e < N; e++) {
        printf("%f + i*%f, ", input[e].real, input[e].imag);
    
    }
    printf("]);\n");

    printf("outputFFT = single([");
    for (int e = 0; e < N; e++) {
        printf("%f + i*%f, ", outputFFT[e].real, outputFFT[e].imag);
    
    }
    printf("]);\n");

    printf("outputDFT = single([");
    for (int e = 0; e < N; e++) {
        printf("%f + i*%f, ", outputDFT[e].real, outputDFT[e].imag);
    
    }
    printf("]);\n");


    printf("MAX SUM DIF %f\n", total_dif);

    return 0;
}
