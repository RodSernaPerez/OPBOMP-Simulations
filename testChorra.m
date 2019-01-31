close all
N=3:2:10;


Pe=0:0.05:0.5;

nTests=100000;
results=zeros(length(N),length(Pe));
for l=1:length(N)
    for i=1:length(Pe)
        for j=1:nTests
            x=(rand(1,N(l))>Pe(i));
            if(length(find(x==1))>length(find(x==0)))
                results(l,i)=results(l,i)+1;
            end
        end

    end
end
results=1-results/nTests;
%%
figure
semilogy(Pe,results);
legend(int2str(N'))