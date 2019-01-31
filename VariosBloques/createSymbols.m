function [symbols,sent_bits,Nbits_position,offset,combinations]=createSymbols(N,nActiveBlocks,D,combinations)
%   CREATESYMBOLS creates symbols in a disperse signal
%   createSymbols(N,nBlocksActivos,D) creates the symbols in a vector of N
%   samples, with nActiveBlocks non-null blocks. D shows the samples that
%   belong to each block. It returns too the number of modulated bits,
%   number of modulated bits in the position, the offset vector and the
%   combinations of the number of blocks.

    [numBlocks,size_blocks]=size(D);
       
    NValoresActivos=size(D,2)*nActiveBlocks; %Number of non-null values

    % The number of bits that are in the position depends on the number of
    % combinations for the blocks
    Nbits_position=floor(log2(nCombinaciones(numBlocks,nActiveBlocks)));
    
    % The number of bits that are in the shape
    Nbits_shape=NValoresActivos;
    
    bits_per_symbol=Nbits_position+Nbits_shape;
    
    % Genera los bits
    sent_bits=round(rand(1,bits_per_symbol));

    % Gets the list of combinations (if it is not already computed)
    if(isempty(combinations))
        combinations=listCombinations(1:numBlocks,nActiveBlocks);
    end
    
    % The bits in the position say which blocks are activated
    indices_seleccionados=combinations(bi2de(sent_bits(1:Nbits_position))+1,:);
    active_blocks=D(indices_seleccionados,:);    
    active_blocks=active_blocks(:);

    % The values are 1 if the bit is 0 and 2 if it is 1
    non_null_values=sent_bits(Nbits_position+1:end)+1;
    
    symbols=zeros(1,N);
    symbols(active_blocks)=non_null_values; %Disperse vector

    % Computes the optimum offset
    offset= size_blocks/(3*NValoresActivos-2*N)*ones(1,N);
    
    
    symbols=symbols+offset;
end

function res=nCombinaciones(N,k)
    % Returns the number of way of placing N elements in groups of K elements
    % without any repetition and no importance to the order
    a=(N-k+1):N;
    b=1:k;

    res=prod(a)/prod(b);
end

function y=listCombinations(list, K)
    %Gives the combination of the elements in list taken by K elements
    if size(list,1)==1
        list=list';
    end
    N=length(list);

    
    y=[];
    if K>1
        for i=1:N
            X=listCombinations(list(i+1:end),K-1);
            Y=zeros(size(X,1),K);
            Y(:,2:end)=X;
            Y(:,1)=ones(1,size(Y,1))*list(i);
            y=[y;Y];
        end
    else
        y=list;
    end
end