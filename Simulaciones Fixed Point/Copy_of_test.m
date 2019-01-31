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

Tsim=2.240e-3; %Tiempo de simbolo
Fini=42e3; %Primera frecuencia usada
B=22e3; %Ancho de banda

Fend=Fini+B; %Ultima frecuencia usada
deltaf=1/Tsim; %Separación entre portadoras

B=Fend-Fini; %Ancho de banda

Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

Fs=length(portadoras)/Tsim; %Frecuencia de muestreo que habria que usar sin CS
Fs=2*portadoras(end)*deltaf;
N=length(portadoras); %Longitud del vector disperso

SNR=0; % En dBs
bits=1;
Fs_CS=0.1*Fs; %Frecuencia de muestreo, mucho inferior a la Nyquist gracias al CS

%% CREACION DE LOS BLOQUES
nBlocksActivos=2;
NValoresActivos=8;
size_blocks=NValoresActivos/nBlocksActivos; %Tamaño de los bloques
saltos=size_blocks; %Diferencia entre el comienzo de un bloque y el del siguiente


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
simbolos=zeros(1,N);

bits_en_posicion=floor(log2(factorial(numBlocks)/(factorial(nBlocksActivos)*factorial(numBlocks-nBlocksActivos))));
bits_en_forma=NValoresActivos;

bits_por_simbolo=bits_en_posicion+bits_en_forma;
bits_enviados=round(rand(1,bits_por_simbolo));

Comb=listarCombinaciones(1:numBlocks,nBlocksActivos);

indices_seleccionados=Comb(bi2de(bits_enviados(1:bits_en_posicion))+1,:);

indices_activos=D(indices_seleccionados,:);    
indices_activos=indices_activos(:);

simbolos_utiles=bits_enviados(bits_en_posicion+1:end)+1;
simbolos(indices_activos)=simbolos_utiles; %Vector disperso

offset= size_blocks/(3*NValoresActivos-2*N)*ones(1,N);
simbolos=simbolos+offset;


%% Creación de la matriz y la señal

PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,N);

%La matriz será ruido gaussiano para las posiciones con portadoras y 0 en
%el resto
s=rng(4);
PREM_MOD_MATRIX(portadoras,:)=(randn(Nportadoras,N)+randn(Nportadoras,N)*1i);
rng('shuffle')

%Tiene que ser simetrica conjugada
PREM_MOD_MATRIX=orth([PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))]);

% Matriz de la idft
F=(dftmtx(size(PREM_MOD_MATRIX,1))/sqrt(size(PREM_MOD_MATRIX,1)))^(-1);

% Llamamos X a la matriz que se va a usar 
X=F*PREM_MOD_MATRIX;
X=real(X); %Son reales, pero para que Matlab lo sepa

% Vamos a forzar que despues del diezmado la matriz sea ortonormal para
% cada bloque
indices=1:Fs/Fs_CS:size(X,1); %Indices de las muestras
for i=1:size(D,1)
    A=X(:,D(i,:));
    M=A(indices,:);
    [~,~,V]=svd(M); %V será la matriz que haga que la matriz A*V sea ortogonal en esas filas 
    A=A*V;
    M=M*V;
    I=diag(sqrt(diag(M'*M)));
    %A=A*I^(-1);
    X(:,D(i,:))=A;
end

PREM_MOD_MATRIX=F^(-1)*X;
x=X*simbolos';
A=X;
% x es el vector que se transmite

%% Esto es para ver como queda el espectro de la señal transmitida

X=F^(-1)*x;

if PLOT==1
    f=Fs/1000*(0:length(x)-1)/length(x)-Fs/2000;
    
    figure
    stem((0:length(x)-1)*1/Fs*1e3,x)
    xlabel('Tiempo [ms]');
    title('Señal Discreta');


    figure
    plot(f,fftshift(abs(X)))
    xlabel('Frecuencia (KHz)')
    title('Espectro discreto');
end

%% PASO A CONTINUO
Fscont=30*Fs; %Será la frecuencia de muestreo que simule el continuo
Ncont=floor(Tsim*Fscont); %Numero de muestras de la señal "continua"
x_cont = interp(x,round(Fscont/Fs)); %La señal "continua"

if PLOT==1
    f=Fscont/1000*(0:length(x_cont)-1)/length(x_cont)-Fscont/2000;
    
    figure
    hold on
    plot((0:length(x_cont)-1)*1/Fscont*1e3,x_cont)
    plot((0:Fscont/Fs:length(x_cont)-1)/Fscont*1e3,x,'x');
    legend('Señal continua','Señal discreta');
    title('Señal continua');
    xlabel('Tiempo [ms]');

    figure
    %plot(linspace(-0.5*Fscont,0.5*Fscont,length(x_cont))/1e3,abs(fftshift(fft(x_cont))));
    plot(f,fftshift(abs(fft(x_cont))))
    xlabel('Frecuencia (KHz)')
    title('Espectro continuo');

end

%% CANAL

% Ruido
power_signal=norm(x_cont)^2;

power_noise=10^(-SNR/10)*power_signal;
noise_limited=BandLimitedNoise(length(x_cont),Fscont,Fini,Fend);
noise_limited=noise_limited/norm(noise_limited)*sqrt(power_noise);


x_cont=x_cont+noise_limited';

if PLOT==1
    f=Fscont/1000*(0:length(x_cont)-1)/length(x_cont)-Fscont/2000;
    
    figure
    plot(f,abs(fftshift(fft(x_cont))));
    xlabel('Frecuencia (KHz)')
    title('Espectro continuo con ruido');

end
%% RECEPTOR NYQUIST
% Es el receptor común que sigue el criterio de Nyquist. Solo es para hacer
% alguna prueba simple.

Fs_rec=Fs; %Frecuencia a la que se va a muestrear
x_rec=x_cont(1:Fscont/Fs_rec:end); %Señal muestreada
X_rec=F^(-1)*x_rec; %Su espectro
if PLOT==1
    f=Fs_rec/1000*(0:length(X_rec)-1)/length(X_rec)-Fs_rec/2000;
    
    figure
    stem(f,fftshift(abs(X_rec)))
    hold on
    plot(f,fftshift(abs(PREM_MOD_MATRIX*simbolos')),'x')
    xlabel('Frecuencia (KHz)')
    title(['Espectro enviado y recuperado con muestrador Nyquist a ',num2str(Fs_rec/1000),' KHz'])
end

%% RECEPTOR CS
% Receptor usando CS


x_CS=changeFs(x_cont,Fscont,Fs_CS);


A=A(indices,:); %Cogemos solo las filas correspondientes a las muestras


x_CS=x_CS-A*offset';
% A es la matriz tal que A*simbolos=x_CS, siendo simbolos disperso

% NOTA: las dos sigueintes lineas son las que estiman los simbolos
i=PrimerBloque(x_CS,D,A,nBlocksActivos); %Busca cual es el bloque en el que se transmite
i=unique(i); %Para que lo ordene
m=D(i,:);
m=m(:);
simbolos_rec=SegundoBloque(x_CS,A,m); %Se estiman los simbolos

%Redondeamos los resultados
simbolos_redondeados=round(simbolos_rec);
simbolos_redondeados(simbolos_redondeados<1.5)=1;
simbolos_redondeados(simbolos_redondeados>1.5)=2;

%Comparamos los simbolos extraidos con los que se enviaron
figure
stem(1:length(simbolos),simbolos-offset);
hold on
plot(m,simbolos_rec,'x');
plot(m,simbolos_redondeados,'s');
legend('Enviados','Recibidos','Redondeados');
title(['Comparacion de lo extraido y lo enviado (SNR=',int2str(SNR),'dBs)']);

bits_en_forma_recibidos=simbolos_redondeados-1;

I=ismember(Comb,i);
I=sum(I,2);
[~,I]=max(I);
bits_en_fase_recibidos=de2bi(I-1,bits_en_posicion);
bits_recibidos=[bits_en_fase_recibidos,bits_en_forma_recibidos'];

%% MEDIDAS
capacidad=length(bits_recibidos)/Tsim;
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

errores=length(find(bits_recibidos'~=bits_enviados'));
fprintf('Ha habido %d errores (un %.2f por ciento)\n',errores,100*errores/length(bits_recibidos));


