function [Matriz]=crearMatriz(Exp)

ColumnasOrtogonales=1;

%% PARAMETROS
Tsim=Exp.Tsim; %Tiempo de simbolo
Fini=Exp.Fini; %Primera frecuencia usada
Fend=Exp.Fend; %Última frecuencia usada

N=Exp.Bloques.N;
size_blocks=Exp.Bloques.size_blocks;

deltaf=1/Tsim; %Separación entre portadoras

B=Fend-Fini; %Ancho de banda

Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

Fs=NFFT/Tsim; %Frecuencia de muestreo que habria que usar sin CS

%% MATRIZ
if(ColumnasOrtogonales==0)
    PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,N);
    PREM_MOD_MATRIX(portadoras,:)=randn(Nportadoras,N)+randn(Nportadoras,N)*1i;
    PREM_MOD_MATRIX=[PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))];

    for i=1:size(PREM_MOD_MATRIX,2)
        PREM_MOD_MATRIX(:,i)=PREM_MOD_MATRIX(:,i)/norm(PREM_MOD_MATRIX(:,i));
    end
    % PREM_MOD_MATRIX es una matriz gaussiana y compleja conjugada para
    % asergurar que el resultado de multiplicarla por el vector sea el espectro
    % de una señal real
end

%% MATRIZ

if (ColumnasOrtogonales==1)
    %Vamos a crear al principio solo un bloque y solo la mitad de las
    %portadoras
    PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,size_blocks);


    %La matriz será ruido gaussiano para las posiciones con portadoras y 0 en
    %el resto
    PREM_MOD_MATRIX(portadoras,:)=randn(Nportadoras,size_blocks)+randn(Nportadoras,size_blocks)*1i;

    %Tiene que ser simetrica conjugada
    PREM_MOD_MATRIX=[PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))];

    %Calculamos los indices de las muestras que nos vamos a quedar
    Fs_CS=Exp.Fs_CS*Fs; %Frecuencia de muestreo que se va a usar
    indices=1:Fs/Fs_CS:size(PREM_MOD_MATRIX,1); %Indices de las muestras
    indices=round(indices);
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

    PREM_MOD_MATRIX=[PREM_MOD_MATRIX,zeros(size(PREM_MOD_MATRIX,1),N-size(PREM_MOD_MATRIX,2))];
end

Matriz=PREM_MOD_MATRIX;

