function x=ValuesEstimation(y,A,indices,int_fp,dec_fp)
%   VALUESESTIMATION: return the non-cero values of the solution of the CS
%   problem.
%   ValuesEstimation(y,A,indices) returns the non-cero values of the
%   solution of the CS problem from the y signal using the matrix A knowing
%   those values are found in the positions in the vector indices

    M=A(:,indices); %Takes the coloumns of the non-zero values
    M=pinv(M);
    M_fp=convertToFixPoint(M,int_fp,dec_fp);
    y_fp=convertToFixPoint(y,int_fp,dec_fp);

    x=M_fp*y_fp; %Computes the aproximation of the non-zero values
    x=double(x);
end