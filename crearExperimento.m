function [Exp]=crearExperimento()
%% Crea una estructura con valores por defecto

Exp.Tsim=3e-9; %Tiempo de simbolo
Exp.Fini=5e9; %Primera frecuencia usada
Exp.Fend=9e9; %Última frecuencia usada
Exp.bits=1; %Numero de bits en cada elemento del vector no nulo
Exp.SNR=0; % En dBs
Exp.Fs_CS=0.05;

deltaf=1/Exp.Tsim; %Separación entre portadoras
NFFT=2*floor(Exp.Fend/deltaf);
Exp.Fs=NFFT/Exp.Tsim;

Exp.Bloques=crearBloques();

Exp.Matriz=crearMatriz(Exp);

Exp.bits_por_simbolo=Exp.bits*Exp.Bloques.size_blocks+floor(log2(Exp.Bloques.numBlocks));


end
