function noise=BandLimitedNoise(N,Fs,Fini,Fend)

    f=Fs*(0:N-1)/N;
    noise=zeros(1,(N-2)/2);
    for i=1:length(noise)
        if(f(i+1)>Fini && f(i+1)<Fend)
            noise(i)=randn(1,1)+randn(1,1)*j;
        end
    end
    noise=[0,0,noise,flip(conj(noise))];
    
    noise=real(ifft(noise));
    
end