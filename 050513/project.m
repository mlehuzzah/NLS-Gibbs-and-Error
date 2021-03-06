%  Fourier Spliting Method + RK3
 
clear all
close all
b=.5;
alph=pi/2;
pause on

xmax = pi ;

% 2 xmax= computational window in the Fourier domain
% number of  grid points


nx =  1024;

dx= 2*xmax/nx;               % grid spacing
x= dx*(-nx/2:nx/2-1);        % Fourier frequencies
x = x' ;

L= dx*nx/2;
dk= 2*pi/(2*L);
k= dk*(-nx/2:nx/2-1);



dt = 0.001;
 
time = 0 ;

T_final = 6;
Nstep = round(T_final/dt);

% initial condition

u = 0*x+1;
ln=nx/10;
u(round(nx/2-ln):round(nx/2+ln)) = b*exp(i*alph);
s='b=.5, alpha=pi/2, L=pi/10';

plot(x,abs(u))
title(s)
axis([-4 4 0 2.5])
pause(.8)
drawnow
pause
     
for j= 1:Nstep-1
    time = time+dt; 
    
    % first linear half step ------------------------
    u = fourier2(u,dt);            
    % nonlinear full step - TVD RK3n -------------------
  
    %RK1
    dflux = 2*i*abs(u).^2 .* u; 
    u1 = u + dt * dflux;

    %RK2
    dflux = 2*i*abs(u1).^2 .* u1; 
    u2 = 3/4*u + 1/4*(u1 + dt * dflux); 
    
    %RK3
    dflux = 2*i*abs(u2).^2 .* u2; 
    u = 1/3*u + 2/3*(u2 + dt * dflux); 
  
    % linear half step ------------------------------
    u  = fourier2(u,dt);  
	
    
    plot(x,abs(u))
    title(s)
	axis([-4 4 0 2.5])
	%pause(.01)
	drawnow
    %mov(j)=getframe(gcf); %gets the whole frame... including axis and title and not just the graph
    %mov(j).cdata = mov(j).cdata(end :-1: 1, :, :); %for some reason my movie is upside down!!!
end


%movie2avi(mov,'wsp1nG.avi');

