%% s=spark(A) Función que recibe una matriz A y devuelve su spark s. 
% OJO: función absurdamente lenta. Si la llamas para una matriz más grande
% de 4x4 el ordenador huele a salchichas.

function s=spark(A)
    [~,N]=size(A);

    finish=0;
    for s=1:N
       C=nchoosek(1:N,s);
       for i=1:length(C(:,1))
           c=C(i,:);
           B=A(:,c);
           r=rank(B);
           if(r<s)
               finish=1;
               break;           
           end
       end
      if(finish==1)
        break;
       end
    end
end