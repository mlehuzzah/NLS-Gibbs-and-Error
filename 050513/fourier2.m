function v= fourier2(u,deltat)

nx= length(u);
dx=.1;

L= dx*nx/2;
dk= 2*pi/(2*L);
k= dk*(-nx/2:nx/2-1);


uhat= fftshift(fft(fftshift(u)));

for ix = 1:nx
        
        kk = k(ix);
        uhat(ix) = exp(-i*kk*kk*deltat/2)* uhat(ix)*filter_exp((ix-nx/2)/nx); 
   
end

v= fftshift(ifft(fftshift(uhat)));
return
