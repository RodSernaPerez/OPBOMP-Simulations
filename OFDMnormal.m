freqs=87:182;

N=length(freqs);
bits=round(rand(1,N));

x=OFDMMod(bits,512,N,0,0);
SNR=10;
SNR_lin=10^SNR/10;
S=norm(x)^2;

NoiseNorm=SNR_lin/S;
n=
