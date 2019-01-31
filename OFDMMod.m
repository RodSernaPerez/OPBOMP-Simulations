function [ofdmSymbolsSe]=OFDMMod(modSymbols,NFFT,Nf,Ncyc,faseinit)

modSymbols = modSymbols(:).'; % lo pongo en fila

if mod(length(modSymbols),Nf)~=0
    error('modSymbols'' length must be an integer number of tines Nf') ;
end

Nsym = length(modSymbols)/Nf; %Número de símbolos OFDM

% Agrupo en una matriz con tantas filas como simbolos y tantas columnas como subportadoras de datos
modSymbolsBlock = reshape(modSymbols,Nf,Nsym)';   

% Construyo el espectro de OFDM. Utilizo la matriz anterior en las portadoras de datos y pongo ceros en las portadoras que no llevan datos
ofdmSymbolsFreq = zeros(Nsym,NFFT/2);  

freqs=87:182; %Rango de portadoras para datos
pilot=86; %Portadora del piloto

%Creacion de la parte positiva del espectro
ofdmSymbolsFreq(:,freqs)=modSymbolsBlock; %Datos
ofdmSymbolsFreq(:,pilot)=1-2*faseinit(1:Nsym); %Piloto

%Lo duplicamos de forma que quede simetrica conjugada
aux=ofdmSymbolsFreq(:,NFFT/2:-1:1);
aux=conj(aux);
aux=[aux(:,end),aux(:,1:end)];
ofdmSymbolsFreq=[ofdmSymbolsFreq,aux];
ofdmSymbolsFreq=ofdmSymbolsFreq(:,1:end-1);

%Añadimos un ultimo simbolo de ceros
ofdmSymbolsFreq=[ofdmSymbolsFreq;zeros(1,NFFT)];

% Utilizo la IFFT para pasar al dominio del tiempo 
ofdmSymbolsPa = (NFFT/sqrt((Nf)*2))*ifft(ofdmSymbolsFreq.');  

%Añadimos los bits ciclicos
ofdmSymbolsPa=[ofdmSymbolsPa(NFFT-Ncyc+1:NFFT,:);ofdmSymbolsPa];



% Hago el paralelo/serie para sacar las cosas al tiemp
ofdmSymbolsSe = ofdmSymbolsPa(:); 
