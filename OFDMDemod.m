function [modSymbols,fases]=OFDMDemod(ofdmSymbolsSe,NFFT,Nf,Ncyc,Equaliz)


if mod(length(ofdmSymbolsSe),NFFT+Ncyc)~=0
    error('ofdmSymbolsSe'' length must be an integer number of tines NFFT') ;
end

Nsym = length(ofdmSymbolsSe)/(NFFT+Ncyc);%Número de simbolos OFDM

y_symbol = reshape(ofdmSymbolsSe,NFFT+Ncyc,Nsym); % Paso el vector de entrada a una matriz con tantas filas como la longitud de la FFT y tantas columnas como simbolos
y_symbol=y_symbol(Ncyc+1:end,:); %Quito las muestras cíclicas añadidas
Y_symbol = fft(y_symbol)*(NFFT/sqrt((Nf)*2));  % Paso a la frecuencia

Y_symbol=Y_symbol(:,1:end-1);%Quitamos el ultimo simbolo de ceros
Nsym=Nsym-1;


if(length(Equaliz)>1)%Ecualizacion de la señal
    for i=1:Nsym
        Y_symbol(:,i)=Equaliz.'.*Y_symbol(:,i);
    end
end

freqs=87:182; %Rango de portadoras para datos
pilot=86; %Portadora del piloto


yMod_symbols = Y_symbol(freqs,:); % Extraigo las portadoras de datos
fases=angle(Y_symbol(pilot,:)); %Extraigo las fases iniciales 

%modSymbols = reshape(yMod_symbols ,Nsym*Nf,1);  % Vuelvo a convertir en un vector para la variable de salida que ira al demodulador digital.
modSymbols=yMod_symbols(:);
