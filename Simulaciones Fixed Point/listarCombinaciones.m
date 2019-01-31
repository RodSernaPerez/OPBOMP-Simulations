function y=listarCombinaciones(lista, K)

    if size(lista,1)==1
        lista=lista';
    end
    N=length(lista);

    
    y=[];
    if K>1
        for i=1:N
            X=listarCombinaciones(lista(i+1:end),K-1);
            Y=zeros(size(X,1),K);
            Y(:,2:end)=X;
            Y(:,1)=ones(1,size(Y,1))*lista(i);
            y=[y;Y];
        end
    else
        y=lista;
    end

        
    
end