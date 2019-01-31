classdef BLOQUES
    
    properties 
        %Los que se pueden introducir a mano
        N
        size_blocks
        saltos
        
        %Los que no se pueden cambiar
        numBlocks
        D
    end
    
    methods
        function bloques=BLOQUES(Exp)
            bloques.size_blocks=20;
            bloques.N=Exp.N; 
            bloques.saltos=12;
            
            bloques=update(bloques);
        end
        
        function Bloques=update(Bloques)
            Bloques.numBlocks=ceil((Bloques.N-Bloques.size_blocks+1)/Bloques.saltos);
            D2=zeros(Bloques.numBlocks,Bloques.size_blocks);
            D2(1,:)=1:Bloques.size_blocks;
            for i=2:size(D2,1)
                D2(i,:)=D2(i-1,:)+Bloques.saltos;
            end

            Bloques.D=D2;
        end
        
        function Bloques=setN(Bloques,N)
            Bloques.N=N;
            Bloques=update(Bloques);
        end
                
       function Bloques=setsize_blocks(Bloques,size_blocks)
           Bloques.size_blocks=size_blocks;
           Bloques=update(Bloques);
       end        
        
       function Bloques=setsaltos(Bloques,saltos)
            Bloques.saltos=saltos;
            Bloques=update(Bloques);
        end
    end
end