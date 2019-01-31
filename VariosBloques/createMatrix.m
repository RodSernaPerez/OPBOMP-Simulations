function X=createMatrix(N,NFFT,seed,carriers,Fs_Ny,Fs_CS,D)
% CREARMATRIX: crea una matriz para el muestreo de señales dispersas
%   X=crearMatrix(N,NFFT,seed,portadoras,Fs,Fs_CS,D) devuelve la matriz de
%   muestreo. N es el número de muestras de la señal, NFFT es el número de
%   portadoras, seed es la semilla usada para la generación de números
%   aleatorios, portadoras es un vector con las portadoras que se usarán 
%   en OFDM,Fs es la frecuencia de muestro según Nyquist, Fs_CS es la
%   frecuencia de muestreo del receptor de CS y D es la matriz con las
%   portadoras de cada uno de los bloques.

%   CREATEMATRIX: creates a matrix for sampling disperse signals
%   X=createMatrix(N,NFFT,seed,carriers,Fs_Ny,Fs_CS,D) returns the sampling
%   matrix. N is the number of samples of the signal, NFFT is the number
%   of samples that enter the FFT, seed is the seed used for the random
%   generator, carriers is the vector of the carriers that are used for
%   data in OFDM, Fs_Ny is the sampling frequency according to Nyquist,
%   Fs_CS is the sampling frequency in the receiver and D is the matrix
%   with the carriers in each one of the blocks.

    NumberCarriers=length(carriers); %Number of carriers that will be used
    
    %The matrix will be random gaussian variables for the positions of the
    %carriers that are used and  in the others.
    PREM_MOD_MATRIX=zeros((NFFT-2)/2+2,N);
   
    rng(seed); %Initiliziates the generator of random numbers
    
    % Places random complex values to the carriers
    PREM_MOD_MATRIX(carriers,:)=(randn(NumberCarriers,N)+randn(NumberCarriers,N)*1i);
    
    rng('shuffle') %Frees the random generator
    
    %The matrix must be complex conjugated
    PREM_MOD_MATRIX=orth([PREM_MOD_MATRIX;flip(conj(PREM_MOD_MATRIX(2:end-1,:)))]);

    % Matrix of the IDFT
    iF=(dftmtx(size(PREM_MOD_MATRIX,1))/sqrt(size(PREM_MOD_MATRIX,1)))^(-1);

    % Puts the two matrices together
    X=iF*PREM_MOD_MATRIX;
    X=real(X); %They are already real, but MATLAB does not know

    %% Orthogonal Blocks
    % Vamos a forzar que despues del diezmado la matriz sea ortogonales para
    % cada bloque
    % We will force that after the down-sampling the matrix is orthogonal
    % for each block
    
    indices=1:Fs_Ny/Fs_CS:size(X,1); %Indices of the samples that are taken in the receiver

    % For each block
    for i=1:size(D,1)
        
        % A is the submatrix of the block i
        A=X(:,D(i,:));
        
        % M has the rows of the samples that are taken
        M=A(indices,:);
 
        % V is the matrix that will make A*V orthogonal for those rows
        [~,~,V]=svd(M);
        
        A=A*V; % Does the transformation
        
        X(:,D(i,:))=A; % It is inserted in the matrix
    end
end
