function [Bloques]=crearBloques

Bloques.N=100;
Bloques.size_blocks=20;
Bloques.saltos=10;

Bloques.numBlocks=ceil((Bloques.N-Bloques.size_blocks+1)/Bloques.saltos);


D=zeros(Bloques.numBlocks,Bloques.size_blocks);
D(1,:)=1:Bloques.size_blocks;
for i=2:size(D,1)
    D(i,:)=D(i-1,:)+Bloques.saltos;
end

Bloques.D=D;
