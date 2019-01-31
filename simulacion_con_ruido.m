%% PRUEBAS SOBRE FRECUENCIA CON RUIDO
clc
clear
Exp=EXPERIMENT;
numTests=100;
SNR=10;

Fs=1./(15:25);


numBits=Exp.bits_por_simbolo;
 
BER=zeros(length(SNR),length(Fs));
for k=1:length(Fs)
    for m=1:length(SNR)
        Exp.SNR=SNR(m);
        
        Exp=Exp.setFs_CS(Fs(k));
    
        numberErrors=0;
        for i=1:numTests
            bits=round(rand(1,numBits));
            Simbolos=crearSimbolos(bits,Exp);

            Resultados=FuncionTest(Exp,Simbolos);
            aux=decimalToBinaryVector(Resultados.simbolos,Exp.bits)';
            aux=aux(:);
            try
                bits_en_posicion=decimalToBinaryVector(Resultados.posicion-1,Simbolos.bits_posicion);
            catch
                bits_en_posicion=zeros(1,Simbolos.bits_posicion);
            end
            bits_rx=[bits_en_posicion,aux'];
            numberErrors=numberErrors+length(find(round(bits)~=round(bits_rx)));
             close all
        end
        BER(m,k)=numberErrors/(numTests*numBits);
    end
    
end

%% Mostrar imagen
figure
F=surf(SNR,Fs*100,BER);
ylabel('Frecuencia de muestreo/Frecuencia Nyquist (%)');
xlabel('SNR');


%% Guardar imagen

title='Resultados/simulación_con_ruido_';
Numeros=Resultados.numeros;
