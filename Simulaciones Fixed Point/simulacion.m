close all
clc
clear

%% PARAMETERS

Tsym=100e-9; % Symbol length (in seconds)
Fini=2e9; % First frequency used (in Hz)
B=1000e6; % Bandwidth (in Hz)

Fend=Fini+B; % Last frequency used
deltaf=1/Tsym; % Distance between carriers

NumberCarriers=floor(B/deltaf); %Number of carriers used
carriers_indices=floor(Fini/deltaf):floor(Fini/deltaf)+NumberCarriers-1; % indices of the used carriers

NFFT=2*floor(Fend/deltaf); %Number of samples in the FFT

Fs_Nyquist=2*carriers_indices(end)*deltaf; %Nyquist sampling frequency

Fscont=10*Fs_Nyquist; % Sampling frequency to simulate the channel

Fs_CS=0.04*Fs_Nyquist; %Sampling frequency in the CS receiver


% int_fp=[8,16,24];
% dec_fp=[3,7,11];
int_fp=[16];
dec_fp=[7];
%% SIMULATION
numberActiveBlocks=[1]; % Number of active blocks
numberActiveSamples=[3]; % Number of non-null samples

SNR=5; % SNRs to simulate (in dBs)
nTests=1000; %Number of tests

leyend=strings(length(numberActiveSamples),1); %To create the leyend in the graph
BERs=zeros(length(numberActiveSamples),length(SNR)); %To save the measured BERs

capacity=zeros(length(numberActiveSamples),1); %To save the capacities

for k=1:length(int_fp)
    %Creates the matrix that define which samples belong to each block
    D=createBlocks(numberActiveSamples,numberActiveBlocks,NumberCarriers);
    D=D(1:end-1,:);
    
    %Creates the matrix used for sampling
    A=createMatrix(NumberCarriers,NFFT,4,carriers_indices,Fs_Nyquist,Fs_CS,D);
    A_fp=convertToFixPoint(A,int_fp(k),dec_fp(k));
    
    % To save the number of errors in each SNR
    errors=zeros(1,length(SNR));
    
    % To count the number of transmitted bits
    numberTransmitedBits=zeros(1,length(SNR));

    for i=1:length(SNR)
        for j=1:nTests

            % Creates the disperse symbols 
            [symbols,sent_bits,number_bits_in_phase,offset,combinations]=...
                createSymbols(NumberCarriers,numberActiveBlocks,D);

            % Gets the transmitted signal
            x=A*symbols';
            

            %Signal is converted to analogic
            x_cont=Digital2Anal(Fs_Nyquist,Fscont,x);

            % The signal is passed through the channel
            x_received=channel(x_cont,SNR(i),Fscont,Fini,Fend);

            % Bits are extracted from the received signal
            received_bits=receiverCS(x_received,Fscont,Fs_CS,Fs_Nyquist,A,offset,D,combinations,int_fp(k),dec_fp(k));

            % Sums the number of errors
            errors(i)=errors(i)+length(find(received_bits~=sent_bits));
            
            % Sums the number of transmitted bits
            numberTransmitedBits(i)=numberTransmitedBits(i)+length(sent_bits);
           
        end
    end

    %Computes the bit error rates
    BERs(k,:)=errors./numberTransmitedBits;

    % Computes the capacities
    capacity(k)=length(sent_bits)/Tsym;
    
    leyend(k)=strcat(...
        int2str(int_fp(k)),' Bits per word, '...
        , int2str(dec_fp(k)) ,' Bits for decimal part');
    
end
%% PLOTS THE GRAPHS
fig=figure;
semilogy(SNR,BERs,'-x');
legend(leyend');
title(strcat('Sampling frequency= ',int2str(Fs_CS/Fs_Nyquist*100),'% of Nyquist'));

descr = {strcat('Symbol length = ',num2str(round(Tsym*1e9,2)),' ns');
         strcat('Bandwidth = ',num2str(round(B*1e-9,2)),' GHz');
         strcat('CS Sampling = ',num2str(round(Fs_CS*1e-6,2)),' MHz')};
     

annotation('textbox',[0.3 0 .3 .3],'String',descr,'FitBoxToText','on');

xlabel('SNR [dBs]');
ylabel('Probability of error');
grid on
