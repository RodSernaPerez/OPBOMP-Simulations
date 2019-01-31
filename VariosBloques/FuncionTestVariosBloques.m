function [errores,bits_por_simbolo]=FuncionTestVariosBloques(nBlocksActivos,NValoresActivos,SNR,Fs_CS)
%% Hace una simulación de:
% o Generación de los datos a transmitir
% o Modulación para que sea una señal dispersa
% o Transmisión de la señal por el medio AWGN
% o Recuperación de la señal
% o Medición del BER


close all



%% DATOS

Tsim=3e-3; %Tiempo de simbolo
Fini=42e3; %Primera frecuencia usada
B=11e3; %Ancho de banda

Fend=Fini+B; %Ultima frecuencia usada
deltaf=1/Tsim; %Separación entre portadoras

B=Fend-Fini; %Ancho de banda

Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

Fs=length(portadoras)/Tsim; %Frecuencia de muestreo que habria que usar sin CS

N=length(portadoras); %Longitud del vector disperso

Fs_CS=Fs_CS*Fs;


%% CREACION DE LOS BLOQUES

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

%% Crear simbolos


simbolos=zeros(1,N);

bits_en_posicion=floor(log2(factorial(numBlocks)/(factorial(nBlocksActivos)*factorial(numBlocks-nBlocksActivos))));

bits_en_forma=NValoresActivos;

bits_por_simbolo=bits_en_posicion+bits_en_forma;
bits_enviados=round(rand(1,bits_por_simbolo));

Comb=listarCombinaciones(1:numBlocks,nBlocksActivos);

indices_seleccionados=Comb(bi2de(bits_enviados(1:bits_en_posicion))+1,:);

indices_activos=D(indices_seleccionados,:);    
indices_activos=indices_activos(:);


%simbolos_utiles=floor(2^bits*rand(1,NValoresActivos))+1; %Numeros no nulos que se transmiten
simbolos_utiles=bits_enviados(bits_en_posicion+1:end)+1;
simbolos(indices_activos)=simbolos_utiles; %Vector disperso

offset= size_blocks/(3*NValoresActivos-2*N)*ones(1,N);
simbolos=simbolos+offset;


%% Creación de la matriz


PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,N);

%La matriz será ruido gaussiano para las posiciones con portadoras y 0 en
%el resto
rng(4);
PREM_MOD_MATRIX(portadoras,:)=(randn(Nportadoras,N)+randn(Nportadoras,N)*1i);
rng('shuffle')
%Tiene que ser simetrica conjugada
PREM_MOD_MATRIX=orth([PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))]);

% Matriz de la idft
F=(dftmtx(size(PREM_MOD_MATRIX,1))/sqrt(size(PREM_MOD_MATRIX,1)))^(-1);

% Llamamos X a la matriz que se va a usar en el redondeo
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


x=X*simbolos';
A=X;

% x es el vector que se transmite

%% PASO A CONTINUO
Fscont=30*Fs; %Será la frecuencia de muestreo que simule el continuo
x_cont = interp(x,Fscont/Fs); %La señal "continua"


%% CANAL

% Ruido
power_signal=norm(x_cont)^2;

power_noise=10^(-SNR/10)*power_signal;


noise_limited=BandLimitedNoise(length(x_cont),Fscont,Fini,Fend);

noise_limited=noise_limited/norm(noise_limited)*sqrt(power_noise);

x_cont=x_cont+noise_limited';
%% RECEPTOR CS
% Receptor usando CS


x_CS=changeFs(x_cont,Fscont,Fs_CS);


A=A(indices,:); %Cogemos solo las filas correspondientes a las muestras


x_CS=x_CS-A*offset';
% A es la matriz tal que A*simbolos=x_CS, siendo simbolos disperso

% NOTA: las dos sigueintes lineas son las que estiman los simbolos
i=PrimerBloque(x_CS,D,normc(A),nBlocksActivos); %Busca cual es el bloque en el que se transmite
i=unique(i);
m=D(i,:);
m=m(:);
simbolos_rec=SegundoBloque(x_CS,A,m); %Se estiman los simbolos

%Redondeamos los resultados
simbolos_redondeados=round(simbolos_rec);
simbolos_redondeados(simbolos_redondeados<1.5)=1;
simbolos_redondeados(simbolos_redondeados>1.5)=2;

% Extraemos los bits
bits_en_forma_recibidos=simbolos_redondeados-1;
I=ismember(Comb,i);
I=sum(I,2);
[~,I]=max(I);

try 
%Existen algunas combinaciones de bloques imposibles, en ese caso entrará en el catch
% y se pondran todos los bits a cero
    
    bits_en_fase_recibidos=de2bi(I-1,bits_en_posicion);
catch
    bits_en_fase_recibidos=zeros(1,bits_en_posicion);
end
    

bits_recibidos=[bits_en_fase_recibidos,bits_en_forma_recibidos'];



%% MEDIDAS
errores=length(find(bits_recibidos~=bits_enviados));
