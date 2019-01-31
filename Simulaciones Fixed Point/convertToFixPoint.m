function [fpvalue]=convertToFixPoint(value,word_length,decimal_part)
%%  CONVERTTOFIXPOINT: returns the fixed point representation 
%%   convertToFixPoint(value,word_length,decimal_part) 
%%  using word_length bits for each number: 1 bit is for the sign and decimal_part
%% is fot the le
    fpvalue=fi(value,1,word_length,decimal_part);
end