% %% DEPENDENCIA DE LA PROBABILIDAD DE ERROR CON EL SNR
clear
Exp=EXPERIMENT;

numTests=100;
SNR=1:10;

numBits=Exp.bits_por_simbolo;
 
BER=zeros(1,length(SNR));
for k=1:length(SNR)

    numberErrors=0;
    
    for i=1:numTests
        Exp.SNR=SNR(k);
        bits=round(rand(1,numBits));
        Simbolos=crearSimbolos(bits,Exp);

        Resultados=FuncionTest(Exp,Simbolos);
        aux=decimalToBinaryVector(Resultados.simbolos-1,Exp.bits)';
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
    BER(k)=numberErrors/(numTests*numBits);
end

figure
plot(SNR,BER);

% %% PRUEBAS SOBRE FRECUENCIA DE MUESTREO EN AUSENCIA DE RUIDO
% clc
% clear
% Exp=crearExperimento();
% numTests=10;
% SNR=100;
% 
% Fs=1./(90:150);
% % Fs=5.5/100;
% 
% numBits=Exp.bits_por_simbolo;
%  
% BER=zeros(1,length(Fs));
% for k=1:length(Fs)
% 
%     numberErrors=0;
%     
%     for i=1:numTests
%         Exp.SNR=SNR;
%         Exp.Fs_CS=Fs(k);
%         Exp.Matriz=crearMatriz(Exp);
%         bits=round(rand(1,numBits));
%         Simbolos=crearSimbolos(bits,Exp);
% 
%         Resultados=FuncionTest(Exp,Simbolos);
%         aux=decimalToBinaryVector(Resultados.simbolos-1,Exp.bits)';
%         aux=aux(:);
%         try
%             bits_en_posicion=decimalToBinaryVector(Resultados.posicion-1,Simbolos.bits_posicion);
%         catch
%             bits_en_posicion=zeros(1,Simbolos.bits_posicion);
%         end
%         bits_rx=[bits_en_posicion,aux'];
%         numberErrors=numberErrors+length(find(round(bits)~=round(bits_rx)));
%          close all
%     end
%     BER(k)=numberErrors/(numTests*numBits);
% end
% 
% figure
% plot(Fs*100,BER);
% xlabel('Frecuencia de muestreo/Frecuencia Nyquist (%)');
% ylabel('BER');
%%
