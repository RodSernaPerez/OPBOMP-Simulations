function x_cont=Digital2Anal(Fs,Fscont,x)
    x_cont = interp(x,floor(Fscont/Fs)); %La señal "continua"
end