%% Hace una simulación de:
% o Generación de los datos a transmitir
% o Modulación para que sea una señal dispersa
% o Transmisión de la señal por el medio AWGN
% o Recuperación de la señal
% o Medición del BER


close all
clc
clear
PLOT=1;

%% DATOS

Tsim=3e-3; %Tiempo de simbolo
Fini=42e3; %Primera frecuencia usada
%Fend=86e3; %Última frecuencia usada
B=44e3; %Ancho de banda

Fend=Fini+B; %Ultima frecuencia usada
deltaf=1/Tsim; %Separación entre portadoras

B=Fend-Fini; %Ancho de banda

Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

Fs=length(portadoras)/Tsim; %Frecuencia de muestreo que habria que usar sin CS

N=length(portadoras); %Longitud del vector disperso

bits=1; %Numero de bits en cada elemento del vector no nulo

SNR=100; % En dBs

%% CREACION DE LOS BLOQUES

size_blocks=10; %Tamaño de los bloques
saltos=10; %Diferencia entre el comienzo de un bloque y el del siguiente

% NOTA: saltos y size_blocks son diferentes porque podria haber solapes
% entre bloques.

numBlocks=ceil((N-size_blocks+1)/saltos);
D=zeros(numBlocks,size_blocks);
D(1,:)=1:size_blocks;
for i=2:size(D,1)
    D(i,:)=D(i-1,:)+saltos;
end

% D guarda los indices de las posiciones no nulas para cada uno de los
% bloques.

%% MODULACION OFDM

k=size(D,2); %Longitud del bloque
PrimerSimb=D(ceil(size(D,1)*rand),1); %Bloque en el que se modula

simbolos_utiles=floor(2^bits*rand(1,k)); %Numeros no nulos que se transmiten

simbolos=zeros(1,N);

simbolos(PrimerSimb:PrimerSimb+k-1)=simbolos_utiles+1; %Vector disperso


%% Creación de la señal

%Vamos a crear al principio solo un bloque y solo la mitad de las
%portadoras
PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,size_blocks);

s=rng(4);
%La matriz será ruido gaussiano para las posiciones con portadoras y 0 en
%el resto
PREM_MOD_MATRIX(portadoras,:)=randn(Nportadoras,size_blocks)+randn(Nportadoras,size_blocks)*1i;
rng('shuffle')


%Tiene que ser simetrica conjugada
PREM_MOD_MATRIX=[PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))];

%Calculamos los indices de las muestras que nos vamos a quedar
Fs_CS=0.05*Fs; %Frecuencia de muestreo que se va a usar
indices=1:Fs/Fs_CS:size(PREM_MOD_MATRIX,1); %Indices de las muestras

% Matriz de la idft
F=conj(dftmtx(size(PREM_MOD_MATRIX,1)))/size(PREM_MOD_MATRIX,1);

% Llamamos X a la matriz que se va a usar en el redondeo
X=F*PREM_MOD_MATRIX;
X=real(X); %Son reales, pero para que Matlab lo sepa

% Vamos a forzar que la matriz formada por las filas de las muestras que
% nos vamos a quedar sea una matriz ortonormal. De esta forma, la matriz de
% la fase de deteccion y la de estimacion será la misma
X=X(indices,:);
[~,~,V]=svd(X); %V será la matriz que haga que la matriz X*V sea ortogonal en esas filas 

%Se la juntamos a PREM_MOD_MATRIX
PREM_MOD_MATRIX=PREM_MOD_MATRIX*V;

%Ahora hacemos que esa matriz de muestreo ademas tenga columnas normales
X=F*PREM_MOD_MATRIX;
X=real(X);
X=X(indices,:);
I=diag(sqrt(diag(X'*X)));
PREM_MOD_MATRIX=PREM_MOD_MATRIX*I^(-1);

A=PREM_MOD_MATRIX;
for i=1:N/size_blocks-1
    aux=A;
    t=ceil(size(A,2)*rand(1,10));
    aux(:,t)=-A(:,t);
    PREM_MOD_MATRIX=[PREM_MOD_MATRIX,aux];
end
PREM_MOD_MATRIX=[PREM_MOD_MATRIX,zeros(size(PREM_MOD_MATRIX,1),N-size(PREM_MOD_MATRIX,2))];

Mod_symbols=PREM_MOD_MATRIX*simbolos';

% Mod_symbols no es un vector disperso ahora

F=conj(dftmtx(size(PREM_MOD_MATRIX,1)))/size(PREM_MOD_MATRIX,1);
x=F*Mod_symbols;
x=real(x);

% x es el vector que se transmite

%% Esto es para ver como queda el espectro de la señal transmitida

X=fft(x);

if PLOT==1
    figure
    stem((0:length(x)-1)*1/Fs*1e9,x)
    xlabel('Tiempo [ns]');
    title('Señal Discreta');
    mean(X(abs(X)>0.1))
   

    figure
    stem(linspace(-0.5*Fs,0.5*Fs,length(x))/1e9,abs(fftshift(X)))

    xlabel('Frecuencia (GHz)')
    title('Espectro discreto');
end

%% PASO A CONTINUO
Fscont=30*Fs; %Será la frecuencia de muestreo que simule el continuo
Ncont=floor(Tsim*Fscont); %Numero de muestras de la señal "continua"
x_cont = interp(x,Fscont/Fs); %La señal "continua"

if PLOT==1
    figure
    plot((0:length(x_cont)-1)*1/Fscont*1e9,x_cont)
    xlabel('Tiempo [ns]');
    title('Señal continua');
    hold on
    plot((0:Fscont/Fs:length(x_cont)-1)/Fscont*1e9,x,'x');
    legend('Señal continua','Señal discreta');

    figure
    plot(linspace(-0.5*Fscont,0.5*Fscont,length(x_cont))/1e9,abs(fftshift(fft(x_cont))));
    xlabel('Frecuencia (GHz)')
    title('Espectro continuo');
end

%% CANAL

% Ruido
power_signal=norm(x_cont)^2;

power_noise=10^(-SNR/10)*power_signal;
noise=randn(size(x_cont,1),size(x_cont,2));

fn = Fscont;
[B_f, A_f] = fir1(10,[Fini Fend]/fn);
noise_limited = filter(B_f, A_f, noise);
noise_limited=noise_limited/norm(noise_limited)*sqrt(power_noise);

x_cont=x_cont+noise_limited;
%% RECEPTOR NYQUIST
% Es el receptor común que sigue el criterio de Nyquist. Solo es para hacer
% alguna prueba simple.

Fs_rec=Fs; %Frecuencia a la que se va a muestrear
x_rec=x_cont(1:Fscont/Fs_rec:end); %Señal muestreada
X_rec=dftmtx(length(x_rec))*x_rec; %Su espectro
%error=norm(Mod_symbols)/norm(X_rec-Mod_symbols)

%% RECEPTOR CS
% Receptor usando CS

%Fs_CS=0.005*Fs; %Frecuencia de muestreo, mucho inferior a la Nyquist gracias al CS
%x_CS=x_cont(1:Fscont/Fs_CS:end); %Señal muestreada
x_CS=changeFs(x_cont,Fscont,Fs_CS);

A=F*PREM_MOD_MATRIX;
A=real(A);

A=A(1:Fs/Fs_CS:end,:); %Cogemos solo las filas correspondientes a las muestras

% A es la matriz tal que A*simbolos=x_CS, siendo simbolos disperso

% NOTA: las dos sigueintes lineas son las que estiman los simbolos
i=PrimerBloque(x_CS,D,normc(A)); %Busca cual es el bloque en el que se transmite

simbolos_rec=SegundoBloque(x_CS,A,D(i,:)); %Se estiman los simbolos

%Redondeamos los resultados
simbolos_redondeados=round(simbolos_rec);

%Comparamos los simbolos extraidos con los que se enviaron
figure
stem(1:length(simbolos),simbolos);
hold on
plot(D(i,:),simbolos_rec,'x');
plot(D(i,:),simbolos_redondeados,'s');
legend('Enviados','Recibidos','Redondeados');
title(['Comparacion de lo extraido y lo enviado (SNR=',int2str(SNR),'dBs)']);

%% MEDIDAS
capacidad=k*bits*log2(numBlocks)/Tsim;
fprintf('La capacidad es %.2f bits/s\n',capacidad);

fprintf('La frecuencia de muestreo es %.0f Kmuestras/s \n',Fs_CS*10^-3);

num_muestras=length(x_CS);
fprintf('El numero de muestras es %d \n',num_muestras);

tasa_muestras_utiles=(size_blocks+log2(numBlocks))/N;
fprintf('La eficiencia es %.2f \n',tasa_muestras_utiles);

% Degradacion frente al limite de Shanon
SNR_sin_bits_extra=2^(capacidad/B)-1;
SNR_con_bits_extra=2^((capacidad/tasa_muestras_utiles)/B)-1;

degradacion=log10(SNR_con_bits_extra/SNR_sin_bits_extra);
fprintf('La degradacion frente a Shanon es %.2f dBs \n',degradacion);

errores=length(find(simbolos_redondeados'~=simbolos(D(i,:))));
fprintf('Ha habido %d errores (un %.2f por ciento)\n',errores,100*errores/length(D(i,:)));

fprintf('La ganancia en capacidad es %f \n',Fs/Fs_CS*tasa_muestras_utiles);
