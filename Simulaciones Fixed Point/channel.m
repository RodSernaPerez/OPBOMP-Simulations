function y_cont=channel(x_cont,SNR,Fscont,Fini,Fend)
%   CHANNEL: simulates a communication channel as a AWGN channel
%   channel(x_cont,SNR,Fscont,Fini,Fend) returns the results of adding
%   white gaussian noise to x_cont. Fscont is the sampling frecuency used
%   to simulate the continous signal, Fini is the beggining frequency used
%   and Fend is the last one. SNR is the Signal-to-Noise ratio of the
%   channel in dBs of power.

    power_signal=norm(x_cont)^2; %Power of the transmitted channel
    power_noise=10^(-SNR/10)*power_signal; %Power of the noise that must be added
    
    % Crease a vector of noise in the desired band
    noise_limited=BandLimitedNoise(length(x_cont),Fscont,Fini,Fend);
    
    % Changes the power of the vector
    noise_limited=noise_limited/norm(noise_limited)*sqrt(power_noise);

    % Adds the noise to the signal
    y_cont=x_cont+noise_limited';
end

function noise=BandLimitedNoise(N,Fs,Fini,Fend)
% This function returns a vector of N positions with gaussian noise between
% the frecuencies Fini and Fend. Fs is the sampling frecuency of the
% signal.
    f=Fs*(0:N-1)/N; %Vector of frequencies
    
    % Creates the noise in the frequency domain for the desired band
    noise=zeros(1,(N-2)/2);
    for i=1:length(noise)
        if(f(i+1)>Fini && f(i+1)<Fend)
            noise(i)=randn(1,1)+randn(1,1)*j;
        end
    end
    
    %Makes the signal complex conjugated
    noise=[0,0,noise,flip(conj(noise))];
    
    %Converts it to the time domain
    noise=real(ifft(noise));
    
end