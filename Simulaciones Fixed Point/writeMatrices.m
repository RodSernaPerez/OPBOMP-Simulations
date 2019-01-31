load('Experiment.mat');

A=Experiment.Ads;
S=Experiment.S;
D=Experiment.D;
x=Experiment.x;
offset=Experiment.offsetElimination;
%offset=Implementation.offset;

Pseudoinverses=[];
fid=fopen('pseudoinversas.txt','w');
for i=1:size(D,1)
    M=pinv(A(:,D(i,:)));
    M=M';
    M=M(:);
    M_fp=convertToFixPoint(M,16,7);
    M=hex(M_fp');
    M=strrep(strcat(M),' ','');
    fprintf(fid,'%s,\n',M);
end
matrixToVerilog(S,16,7,'MatrixS.txt')
matrixToVerilog(Pseudoinverses,16,7,'MatricesP.txt')

matrixToVerilog(offset,16,7,'offset.txt')

matrixToVerilog(x,16,7,'signal.txt')
