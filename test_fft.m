% Convertir la entrada a single y realizar la FFT
mat_fft = single(fft(Input));

% Calcular las magnitudes de las salidas FFT y DFT
mag_outputFFT = abs(outputFFT);
mag_outputDFT = abs(outputDFT);
mag_mat_fft = abs(mat_fft);

% Calcular las diferencias absolutas y sus sumas
sum_dif_dft_fft = sum(abs(mag_outputFFT - mag_outputDFT));
dif_mat_fft = mag_mat_fft - mag_outputFFT;
dif_mat_dft = mag_mat_fft - mag_outputDFT;
sum_dif_mat_fft = sum(abs(dif_mat_fft));
sum_dif_mat_dft = sum(abs(dif_mat_dft));

% Mostrar los resultados de las sumas de diferencias
disp(['Suma de diferencias (DFT-FFT): ', num2str(sum_dif_dft_fft)]);
disp(['Suma de diferencias (MAT-FFT): ', num2str(sum_dif_mat_fft)]);
disp(['Suma de diferencias (MAT-DFT): ', num2str(sum_dif_mat_dft)]);

% Gráfico de diferencias entre mat_fft y outputDFT
figure;
plot(dif_mat_dft);
title('Diferencias entre magnitudes de mat\_fft y outputDFT');
xlabel('Índice');
ylabel('Diferencia de magnitudes');

% Gráfico de diferencias entre mat_fft y outputFFT
figure;
plot(dif_mat_fft);
title('Diferencias entre magnitudes de mat\_fft y outputFFT');
xlabel('Índice');
ylabel('Diferencia de magnitudes');

% Gráfico de la magnitud de la entrada
figure;
plot(abs(Input));
title('Magnitud de la entrada');
xlabel('Índice');
ylabel('Magnitud');

% Gráfico comparativo de magnitudes
figure;
hold on;
plot(abs(mat_fft), 'DisplayName', 'mat\_fft');
plot(abs(outputFFT), 'DisplayName', 'outputFFT');
plot(abs(outputDFT), 'DisplayName', 'outputDFT');
hold off;
title('Comparación de magnitudes');
xlabel('Índice');
ylabel('Magnitud');
legend;