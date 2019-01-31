function i=PrimerBloque(y,D,A)

    Matrix=zeros(size(D,1),size(A,2));
    for i=1:size(Matrix,1)
        Matrix(i,D(i,:))=1;
    end
    
    S=Matrix*A';
    k=S*y;

    figure
    stem(k)
    title('Correlaciones con los bloques');
    
    [~,i]=max(k);
    
end