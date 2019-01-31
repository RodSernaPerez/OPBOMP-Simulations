function x=ValuesEstimation(y,A,indices)
%   VALUESESTIMATION: return the non-cero values of the solution of the CS
%   problem.
%   ValuesEstimation(y,A,indices) returns the non-cero values of the
%   solution of the CS problem from the y signal using the matrix A knowing
%   those values are found in the positions in the vector indices

    M=A(:,indices); %Takes the coloumns of the non-zero values

    x=pinv(M)*y; %Computes the aproximation of the non-zero values
end