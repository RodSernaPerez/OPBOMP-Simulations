function Simbolos=crearSimbolos(bits,Exp)

    numBlocks=Exp.Bloques.numBlocks;
    bits_por_simbolo=Exp.bits;
    

    Simbolos.bits_posicion=floor(log2(numBlocks)); %Numero de bits que se codifican en la posicion
    
    Simbolos.numBlock=bi2de(bits(1:Simbolos.bits_posicion),'left-msb')+1;
    bits=bits(Simbolos.bits_posicion+1:end);
    Simbolos.numeros=bi2de(reshape(bits,length(bits)/bits_por_simbolo,bits_por_simbolo),'left-msb')+1;

end