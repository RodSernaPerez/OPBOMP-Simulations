function matrixToVerilog(A,bits,bits_decimal,namefile)

    [M,N]=size(A);

    total_bits=M*N*bits;

    fid=fopen(namefile,'w');
    k=total_bits;
    for i=1:M
        a=A(i,:);
        a_fp=convertToFixPoint(a,bits,bits_decimal);
        a_fp=a_fp(:)';
        a_fp=hex(a_fp);
        for j=1:length(a_fp)
            a_fp=char(a_fp);
        end
        fprintf(fid,'assign out[%d:%d]=%d''h%s;\n',k-1,k-bits*N,bits*N,strrep(strcat(a_fp),' ',''));
        k=k-bits*N;
    end
    fclose(fid);
end