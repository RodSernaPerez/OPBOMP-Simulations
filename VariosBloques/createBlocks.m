function D=createBlocks(NumberActiveSamples,NumberActiveBlocks,N)
% CREARBLOQUES Devuelve una matriz donde la fila i contiene las portadoras asociadas al
% bloque i.
%   D=CREARBLOQUES(NValoresActivos,nBlocksActivos,N) devuelve una matriz
%   donde la fila i contiene las portadoras asociadas al bloque i.
%   NValoresActivos es el número de portadoras activas, nBlocksActivos es el
%   número de bloques que puede haber activos a la vez y N es el número de
%   portadoras.

%   CREATEBLOCKS Returns the matrix where the row i has the indices of the
%   carriers asociated to the block i. 
%   createBlocks(NumberActiveSamples,NumberActiveBlocks,N) Returns the 
%   matrix where the row i has the indices of the carriers asociated to the
%   block i. NumberActiveSamples is the total number of non-zero samples,
%   NumberActiveBlocks is the number of not null blocks and N is the length
%   of the signal.

    size_blocks=NumberActiveSamples/NumberActiveBlocks; %Size of the blocks
    jumps=size_blocks; %Difference between one block and the next one
    
    % NOTE: jumps and size_blocks can be different if there are
    % intersection between different groups.

    numBlocks=ceil((N-size_blocks+1)/jumps); %Number of different blocks
    
    %Creates the matrix
    D=zeros(numBlocks,size_blocks);
    D(1,:)=1:size_blocks;
    for i=2:size(D,1)
        D(i,:)=D(i-1,:)+jumps;
    end
end