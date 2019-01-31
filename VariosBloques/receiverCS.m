function recived_bits=receiverCS(x,Fs_x,Fs_CS,Fs,A,offset,D,Comb)
%   RECEIVERCS: returns the decoded bits using a Compressive Sampling decoder
%   receiverCS(x_cont,Fscont,Fs_CS,Fs,A,offset,D,numberActiveBlocks,Comb,numberBitsInPosition)
%   returns the bits decoded with the Compressive Sampling demodulator from
%   the signal x, sampled with frequency Fs_x. Fs_CS is the sampling
%   frequency in the receiver, Fs is the sampling frequency used in the
%   transmitter, A is the matrix used to modulate the sparse signal, offset
%   is the vector that is sum to the disperse signal before being modulated
%   D is the matrix which contains the indices of the samples that belong
%   to each block, Comb is the list of combinations of the blocks
%   that is used to modulate in the phase.
   
numberActiveBlocks=size(Comb,2); %Number of blocks that can be actived 
numberBitsInPosition= floor(log2(size(Comb,1))); %Number of bits modulated in the phase

indices=1:Fs/Fs_CS:size(A,1); %Indices of the rows of A that corresponde to the picked samples

x_CS=changeFs(x,Fs_x,Fs_CS); %Downsamples the vector x to the frequency of the CS receiver

A=A(indices,:); % Takes the rows of the matrix of the picked samples

x_CS=x_CS-A*offset'; % Undo the effects of the offset

% Selects the indeces of the blocks that are not null
indices_of_blocks=DetectionModule(x_CS,D,normc(A),numberActiveBlocks); 
indices_of_blocks=unique(indices_of_blocks); %Sorts the indices

%Gets the indices of the samples that correspond to the non-null samples
indices_of_samples=D(indices_of_blocks,:); 
indices_of_samples=indices_of_samples(:);

%Symbols are estimeted
recived_symbols=ValuesEstimation(x_CS,A,indices_of_samples);

% Symbols are rounded to the closest possiblity
round_symbols=round(recived_symbols);
round_symbols(round_symbols<1.5)=1;
round_symbols(round_symbols>1.5)=2;

% Bits from the form are extracted
bits_in_form=round_symbols-1;

% From the indices of the blocks, the bits are extracted using the list of
% combinations
indices_of_blocks=ismember(Comb,indices_of_blocks);
indices_of_blocks=sum(indices_of_blocks,2);
[~,indices_of_blocks]=max(indices_of_blocks);

try 
    % Bits are exptracted from the phase
    bits_in_phase=de2bi(indices_of_blocks-1,numberBitsInPosition);
catch
    % There are some combinations of blocks that are not possible, in this
    % case the bits are set to 0
    bits_in_phase=zeros(1,numberBitsInPosition);
end

recived_bits=[bits_in_phase,bits_in_form'];
end