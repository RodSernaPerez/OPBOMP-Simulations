length_block=3;
number_blocks=8;

Total_points=length_block*number_blocks;

values=0:2^(log2(number_blocks))-1;

values_bits=de2bi(values);
figure

subplot(length(values)/2,2,1);
i=1;

for j=1:length(values)
    x=zeros(1,Total_points);
    
    startingPoint=bi2de(values_bits(j,1:log2(number_blocks)))*length_block+1;
    valores=[1,1,1];
    x(startingPoint:startingPoint+length_block-1)=valores;
    subplot(length(values)/2,2,j);
    
    stem(x);
    ylim([-1.5,1.5])
    title(['\color{red}',mat2str(values_bits(j,:)),'\color{black}',mat2str(valores)]);
   
end

%%
length_block=3;
number_blocks=8;

Total_points=length_block*number_blocks;

values=0:2^(length_block)-1;

values_bits=de2bi(values);
figure

subplot(length(values)/2,2,1);
i=1;

for j=1:length(values)
    x=zeros(1,Total_points);
    
    startingPoint=3*length_block-1;
    valores=values_bits(j,:)*2-1;
    
    x(startingPoint:startingPoint+length_block-1)=valores;
    subplot(length(values)/2,2,j);
    
    stem(x);
    ylim([-1.5,1.5])
    title(['\color{red}',mat2str([0,1,0]),'\color{black}',mat2str(values_bits(j,:))]);
   
end