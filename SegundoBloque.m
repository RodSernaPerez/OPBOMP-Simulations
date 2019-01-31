function x=SegundoBloque(y,A,indices)

    M=A(:,indices);

    x=pinv(M)*y;
end