function [Resultados]=FuncionTest(Exp,Simbolos)

    % persistent F
    F=[];

    if nargin==0
        Exp=EXPERIMENT;
        Exp.SNR=10;
        Simbolos=crearSimbolos(round(rand(1,Exp.bits_por_simbolo)),Exp);
    end

    PLOT=1;
    TEX=0;

    %% PARAMETROS
    Tsim=Exp.Tsim; %Tiempo de simbolo
    Fini=Exp.Fini; %Primera frecuencia usada
    Fend=Exp.Fend; %Última frecuencia usada


    bits=Exp.bits; %Numero de bits en cada elemento del vector no nulo

    SNR=Exp.SNR; % En dBs


    deltaf=1/Tsim; %Separación entre portadoras

    B=Fend-Fini; %Ancho de banda

    Nportadoras=floor(B/deltaf); %Numero de portadoras que se usan

    portadoras=floor(Fini/deltaf):floor(Fini/deltaf)+Nportadoras-1; %indices de las portadoras

    NFFT=2*floor(Fend/deltaf); %Numero de puntos de la FFT (distinto de Nportadoras porque no se usan todas las portadoras)

    Fs=NFFT/Tsim; %Frecuencia de muestreo que habria que usar sin CS


    %% BLOQUES

    N=Exp.Bloques.N; %Longitud del vector disperso
    size_blocks=Exp.Bloques.size_blocks; %Tamaño de los bloques
    saltos=Exp.Bloques.saltos; %Diferencia entre el comienzo de un bloque y el del siguiente

    % NOTA: saltos y size_blocks son diferentes porque podria haber solapes
    % entre bloques.

    D=Exp.Bloques.D;

    numBlocks=Exp.Bloques.numBlocks;
    %% MODULACION OFDM

    k=size(D,2); %Longitud del bloque

    PrimerSimb=D(Simbolos.numBlock,1); %Bloque en el que se modula

    simbolos_utiles=Simbolos.numeros-1; %Numeros no nulos que se transmiten

    simbolos=zeros(1,N);

    simbolos(PrimerSimb:PrimerSimb+k-1)=simbolos_utiles+1; %Vector disperso
    
    offset= size_blocks/(3*size_blocks-2*N)*ones(1,N);
    simbolos=simbolos+offset;

    %% MATRIZ
    PREM_MOD_MATRIX=Exp.Matriz.Matriz;

    Mod_symbols=PREM_MOD_MATRIX*simbolos';

    if(isempty(F))%Asi solo la calcula una vez
        F=sqrt(length(Mod_symbols))*dftmtx(NFFT)^(-1); %Matriz de la ifft
    end

    x=F*Mod_symbols;
    x=real(x);

    %% PASO A CONTINUO
    % Fscont=30*Fs; %Será la frecuencia de muestreo que simule el continuo
    % Ncont=floor(Tsim*Fscont); %Numero de muestras de la señal "continua"
    % x_cont = interp(x,Fscont/Fs); %La señal "continua"
    x_cont=x;
    Fscont=Fs;
    %% CANAL

    % Ruido
    power_signal=norm(x)^2;

    power_noise=10^(-SNR/10)*power_signal;
    noise=randn(size(x,1),size(x,2));

    fn = Fscont;
    [B_f, A_f] = fir1(10,[Fini Fend]/fn);
    noise_limited = filter(B_f, A_f, noise);
    noise_limited=noise_limited/norm(noise_limited)*sqrt(power_noise);

    x=x+noise_limited;


    %% RECEPTOR CS
    % Receptor usando CS


    Fs_CS=Exp.Fs_CS*Fs;
    %x_CS=x_cont(1:round(Fscont/Fs_CS):end); %Señal muestreada

    x_CS=changeFs(x,Fs,Fs_CS);

    A=F*PREM_MOD_MATRIX;
    A=real(A);
    A=A(1:Fs/Fs_CS:end,:);

    % A es la matriz tal que A*simbolos=x_CS, siendo simbolos disperso

    %Quitamos el offset
    x_CS=x_CS-A*offset';
    
    % NOTA: las dos sigueintes lineas son las que estiman los simbolos
    i=PrimerBloque(x_CS,D,normc(A)); %Busca cual es el bloque en el que se transmite
    simbolos_rec=SegundoBloque(x_CS,A,D(i,:)); %Se estiman los simbolos

    %Redondeamos los resultados
    simbolos_redondeados=round(simbolos_rec);
    simbolos_redondeados(simbolos_redondeados<1.5)=1;
    simbolos_redondeados(simbolos_redondeados>2.5)=2;

    if (PLOT==1)
        %Comparamos los simbolos extraidos con los que se enviaron
        figure
        stem(1:length(simbolos),simbolos-offset);
        hold on
        plot(D(i,:),simbolos_rec,'x');
        plot(D(i,:),simbolos_redondeados,'s');
        legend('Enviados','Recibidos','Redondeados');
        title(['Comparacion de lo extraido y lo enviado (SNR=',int2str(SNR),'dBs)']);
    end

    %% MEDIDAS
    if(TEX==1)
        capacidad=k*bits*log2(numBlocks)/Tsim;
        fprintf('La capacidad es %.2f Gbits/s\n',capacidad*10^-9);

        fprintf('La frecuencia de muestreo es %.0f kmuestras/s \n',Fs_CS*10^-6);

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

    end

    Resultados.posicion=i;
    Resultados.simbolos=simbolos_redondeados;
    Resultados.capacidad=k*bits*log2(numBlocks)/Tsim;
end