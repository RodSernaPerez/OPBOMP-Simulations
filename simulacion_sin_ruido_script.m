%% PRUEBAS SOBRE FRECUENCIA DE MUESTREO EN AUSENCIA DE RUIDO
% Dibuja la capacidad obtenida en el enlace para distintos muestreadores.
% La gráfica que se dibuja corresponde a Capacidad*(0.5-Perror)

clc
clear


Exp=EXPERIMENT;
numTests=100;
Exp.SNR=10000;

Fs=1./(5:35);


numBits=Exp.bits_por_simbolo;
 
BER=zeros(1,length(Fs));
for k=1:length(Fs)

    numberErrors=0;
    Exp=Exp.setFs_CS(Fs(k));
    
    for i=1:numTests
        bits=round(rand(1,numBits));
        Simbolos=crearSimbolos(bits,Exp);

        Resultados=FuncionTest(Exp,Simbolos);
        aux=decimalToBinaryVector(Resultados.simbolos,Exp.bits)';
        aux=aux(:);
        try
            bits_en_posicion=decimalToBinaryVector(Resultados.posicion-1,Simbolos.bits_posicion);
        catch
            bits_en_posicion=zeros(1,Simbolos.bits_posicion);
        end
        bits_rx=[bits_en_posicion,aux'];
        numberErrors=numberErrors+length(find(round(bits)~=round(bits_rx)));
         close all
    end
    BER(k)=numberErrors/(numTests*numBits);
end


%% Dibujo de las gráficas

figure
hold on;
xlabel('Frecuencia de muestreo/Frecuencia Nyquist (%)');
ylabel('BER');
title(['SNR=',int2str(Exp.SNR)]);
MinimoFs=Exp.Bloques.size_blocks/(Exp.Fs*Exp.Tsim)*100;
semilogy(Fs*100,BER); %Dibujo de la probabilidad de error con CS
set(gca, 'YScale', 'log')

SNR_lin=10^(Exp.SNR/10); % SNR en escala lineal

% La capacidad es el ancho de banda multiplicado por el nbits por
% portadora, en este caso 2 (2*Fs*Exp.Fs). Al multiplicar por el Tsim
% sacamos el numero de bis
Bits_OFDM_QPSK=2*Fs*Exp.Fs*Exp.Tsim; % Capacidad de QPSK
Bits_OFDM_QPSK=Fs*Exp.Fs*Exp.Tsim/2; % Capacidad de QPSK

SNR_bit=SNR_lin./Bits_OFDM_QPSK; %SNRs por bit transmitido (vector)
SNR_bit=SNR_lin.*Bits_OFDM_QPSK; %SNRs por bit transmitido (vector)

BER_QPSK=qfunc(sqrt(2*SNR_bit)); %Probabilidad de error usando QPSK


Bits_OFDM_BPSK=Fs*Exp.Fs*Exp.Tsim; %Numero de bits en BPSK tansmistidos
SNR_bit=SNR_lin./Bits_OFDM_BPSK; %SNRs por bit transmitido (vector)

BER_BPSK=qfunc(sqrt(2*SNR_bit)); %Probabilidad de error usando BPSK
semilogy(Fs*100,BER_BPSK);
semilogy(Fs*100,BER_QPSK);
legend('Con CS','BPSK','QPSK');
%%
figure
hold on 
title(['Capacidad con SNR=',int2str(Exp.SNR)]);
xlabel('Frecuencia de muestreo (KHz)');
ylabel('KBits/s');

plot(Fs*Exp.Fs/1000,Resultados.capacidad*(1-BER)/1000); %Dibujo para CS
plot(Fs*Exp.Fs/1000,Fs*Exp.Fs.*(1-BER_BPSK)/1000); %Dibujo para BPSK
plot(Fs*Exp.Fs/1000,2*Fs*Exp.Fs.*(1-BER_QPSK)/1000); %Dibujo para QPSK
legend('CS','BPSK','QPSK');
grid on





