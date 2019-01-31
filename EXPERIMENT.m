classdef EXPERIMENT
    
    properties 
        %Los que se pueden introducir a mano
        Tsim %Tiempo de simbolo
        Fini %Primera frecuencia usada
        B %Última frecuencia usada
        bits %Numero de bits en cada elemento del vector no nulo
        SNR% En dBs
        Fs_CS
        N
        
        
        %Los que no
        deltaf %Separación entre portadoras
        NFFT
        Fend
        Fs
        Bloques
        Matriz
        bits_por_simbolo
        
    end
    
    methods
        function Exp=EXPERIMENT()
            Exp.Tsim=3e-3; %Tiempo de simbolo
            Exp.Fini=42e3; %Primera frecuencia usada
            Exp.B=44e3;
            Exp.bits=1; %Numero de bits en cada elemento del vector no nulo
            Exp.SNR=100; % En dBs
            Exp.Fs_CS=0.05;
            Exp.N=100;
            
            Exp=update(Exp);
        end
        
        function Exp=update(Exp)
            Exp.deltaf=1/Exp.Tsim; %Separación entre portadoras
            Exp.NFFT=2*floor(Exp.Fend/Exp.deltaf);
            Exp.Fs=Exp.NFFT/Exp.Tsim;
            Exp.Fend=Exp.Fini+Exp.B;
                   
            Exp.Bloques=BLOQUES(Exp);
            Exp.bits_por_simbolo=Exp.bits*Exp.Bloques.size_blocks+floor(log2(Exp.Bloques.numBlocks));
            Exp.Matriz=MATRIX(Exp);
        end  
    
        function Exp=setTsim(Exp,time)
            Exp.Tsim=time;
            Exp=update(Exp);
        end
        
        function Exp=setFini(Exp,Fini)
            Exp.Fini=Fini;
            Exp=update(Exp);
        end
        
        function Exp=setFend(Exp,Fend)
            Exp.Fend=Fend;
            Exp=update(Exp);
        end
        
        function Exp=setbits(Exp,bits)
            Exp.bits=bits;
            Exp=update(Exp);
        end
        
        function Exp=setSNR(Exp,SNR)
            Exp.SNR=SNR;
            Exp=update(Exp);
        end
        
        function Exp=setBloques(Exp,Bloques)
            Exp.Bloques=Bloques;
            Exp=update(Exp);
        end
        
        function Exp=setFs_CS(Exp,Fs_CS)
            Exp.Fs_CS=Fs_CS;
            Exp=update(Exp);
        end
    end
end
        