function res=nCombinaciones(N,k)
%% Calcula el numero de formas de colocar N elementos en grupos de K elementos,
%% sin repetir ninguno y sin importar el orden.
a=(N-k+1):N;
b=1:k;

res=prod(a)/prod(b);
end