function testFrecuenciaMuestreo(p)
%% DATOS

Tsim=50e-9; %Tiempo de simbolo
Fini=5e9; %Primera frecuencia usada
Fend=9e9; %Última frecuencia usada

deltaf=1/Tsim; %Separación entre portadoras

B=Fend-Fini; %Ancho de banda

Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

Fs=length(portadoras)/Tsim; %Frecuencia de muestreo que habria que usar sin CS

N=length(portadoras); %Longitud del vector disperso

bits=1; %Numero de bits en cada elemento del vector no nulo

SNR=1000; % En dBs

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

simbolos_utiles=floor(2^bits*rand(1,k))+1; %Numeros no nulos que se transmiten

simbolos=zeros(1,N);

simbolos(PrimerSimb:PrimerSimb+k-1)=simbolos_utiles; %Vector disperso


%% Creación de la señal
% 
% PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,N);
% PREM_MOD_MATRIX(portadoras,:)=randn(Nportadoras,N)+randn(Nportadoras,N)*1i;
% PREM_MOD_MATRIX=[PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))];
% 
% for i=1:size(PREM_MOD_MATRIX,2)
%     PREM_MOD_MATRIX(:,i)=PREM_MOD_MATRIX(:,i)/norm(PREM_MOD_MATRIX(:,i));
% end
% PREM_MOD_MATRIX es una matriz gaussiana y compleja conjugada para
% asergurar que el resultado de multiplicarla por el vector sea el espectro
% de una señal real


% OTRA FORMA DE HACER LA MATRIZ MAS EFICIENTE

%Vamos a crear al principio solo un bloque y solo la mitad de las
%portadoras
PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,size_blocks);


%La matriz será ruido gaussiano para las posiciones con portadoras y 0 en
%el resto
PREM_MOD_MATRIX(portadoras,:)=randn(Nportadoras,size_blocks)+randn(Nportadoras,size_blocks)*1i;

%Tiene que ser simetrica conjugada
PREM_MOD_MATRIX=[PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))];

%Calculamos los indices de las muestras que nos vamos a quedar
Fs_CS=p*Fs; %Frecuencia de muestreo que se va a usar
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

%Ahora hacemos que esa matriz de muestreo ademas tenga columnas ortogonales
X=F*PREM_MOD_MATRIX;
X=real(X);
X=X(indices,:);
I=diag(sqrt(diag(X'*X)));
PREM_MOD_MATRIX=PREM_MOD_MATRIX*I^(-1);

A=PREM_MOD_MATRIX;
% t=size(A,2)/(N/size_blocks);
for i=1:N/size_blocks-1
    aux=A;
%     aux(:,t*(i-1)+1:t*(i-1)+t)=-A(:,t*(i-1)+1:t*(i-1)+t);
    t=ceil(size(A,2)*rand(1,10));
    aux(:,t)=-A(:,t);
    PREM_MOD_MATRIX=[PREM_MOD_MATRIX,aux];
end
%PREM_MOD_MATRIX=eye(size(PREM_MOD_MATRIX));
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

%Filtro
% T_max_filtro=0.2*Tsim;
% N_muestras_filtro=round(T_max_filtro*Fscont);
% 
% h=zeros(1,length(x_cont));
% h(1:N_muestras_filtro)=rand(1,N_muestras_filtro)-0.5;
% 
% x_cont_filtrado=conv(x_cont,h);
% x_cont=x_cont_filtrado(1:length(x_cont))';

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
x_CS=x_cont(1:Fscont/Fs_CS:end); %Señal muestreada


% h2=h(1:Fscont/Fs_rec:end);
% H=zeros(length(h2));
% for i=size(H,1):-1:1
%     H(i,:)=h2;
%     h2=[h2(2:end),0];
% end
A=F*PREM_MOD_MATRIX;
A=real(A);

A=A(1:Fs/Fs_CS:end,:); %Cogemos solo las filas correspondientes a las muestras

% A es la matriz tal que A*simbolos=x_CS, siendo simbolos disperso

% NOTA: las dos sigueintes lineas son las que estiman los simbolos
i=PrimerBloque(x_CS,D,normc(A)); %Busca cual es el bloque en el que se transmite

simbolos_rec=SegundoBloque(x_CS,A,D(i,:)); %Se estiman los simbolos

%Redondeamos los resultados
simbolos_redondeados=round(simbolos_rec);
simbolos_redondeados(simbolos_redondeados<1)=1;
simbolos_redondeados(simbolos_redondeados>2^bits)=2^bits;

%Comparamos los simbolos extraidos con los que se enviaron
figure
stem(1:length(simbolos),simbolos);
hold on
plot(D(i,:),simbolos_rec,'x');
plot(D(i,:),simbolos_redondeados,'s');
legend('Enviados','Recibidos','Redondeados');
title(['Comparacion de lo extraido y lo enviado (SNR=',int2str(SNR),'dBs)']);

%% MEDIDAS
clc
capacidad=k*bits*log2(numBlocks)/Tsim;
fprintf('La capacidad es %.2f Gbits/s\n',capacidad*10^-9);

fprintf('La frecuencia de muestreo es %.0f Mmuestras/s \n',Fs_CS*10^-6);

num_muestras=length(x_CS);
fprintf('El numero de muestras es %d \n',num_muestras);

tasa_muestras_utiles=size_blocks/N;
fprintf('La eficiencia es %.2f \n',tasa_muestras_utiles);

% Degradacion frente al limite de Shanon
SNR_sin_bits_extra=2^(capacidad/B)-1;
SNR_con_bits_extra=2^((capacidad/tasa_muestras_utiles)/B)-1;

degradacion=log10(SNR_con_bits_extra/SNR_sin_bits_extra);
fprintf('La degradacion frente a Shanon es %.2f dBs \n',degradacion);

errores=length(find(simbolos_redondeados'~=simbolos(D(i,:))));
fprintf('Ha habido %d errores (un %.2f por ciento)\n',errores,100*errores/length(D(i,:)));

const=0.1;
SNR_noCS=2^(1/const)-1;
SNR_CS=2^(tasa_muestras_utiles*capacidad/(const*B))-1;
SNR_CS/SNR_noCS