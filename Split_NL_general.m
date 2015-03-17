%
%
%  Fourier Spliting Method + RK3
%
%

 
clear all
close all

inx_L = 0;    % Domain index 

for iL = 1:1

inx_L = inx_L + 1;   %interval [0, 1] 
    
xmax = pi; 

% 2 xmax= computational window in the Fourier domain
% number of  grid points


nx =  1024;

dx= 2*xmax/nx;               % grid spacing
x= dx*(-nx/2:nx/2-1);        % Fourier frequencies
x = x'; 

L= dx*nx/2;
dk= 2*pi/(2*L);
k= dk*(-nx/2:nx/2-1);
n0 = nx/2 + 1; 

i = sqrt(-1);

dt = 0.0005;
 
time = 0; 

T_final = 400;
Nstep = round(T_final/dt); 

% initial condition - need to be changed 
u = 0*x + 1;

u(nx/4:nx/2) = 2; 

%----------------------------------
% Derivative matrix

ddx = 2*pi/nx; xx = zeros(nx,1); 
for ix = 1:nx
    xx(ix) = 0 + (ix-1)*ddx; 
end
    
D = zeros(nx);
D2 = zeros(nx); 

for ix = 1:nx
    xt = xx(ix); 
    for iy = 1:nx
        yt = xx(iy); 
        arg = (xt-yt)/2; 
        
        if ix == iy
            D(ix,iy) = 0; 
            D2(ix,iy) = -(nx^2+2)/12;
        else  
            D(ix,iy) = 0.5*cos(nx*arg)*cot(arg) ... 
                      -0.5/nx*sin(nx*arg)/((sin(arg))^2);
            D2(ix,iy) = -0.5*(-1)^(ix+iy)/((sin(arg))^2);
        end
    end
end

D = 2*pi/(xmax*2)*D; 
D2 = (2*pi/(xmax*2))^2 * D2; 

% Interpolating polynomials

Np = 1000;
xi = linspace(0,2*pi,Np);
xp = xi/2/pi*(2*xmax)-xmax; 
xp = linspace(-xmax,xmax,Np); 
xc = (x + xmax)/(2*xmax)*2*pi; 
xc = linspace(0,2*pi,nx+1); 
g = zeros(nx,Np); 
for ix = 1:nx
    g(ix,:) = sin(nx*(xi - xc(ix))/2).*cot((xi - xc(ix))/2)/nx; 
end


% First linear half step


alpha = 100; 
ep = 0.036; 
ep2 = ep*ep;

% first linear half step ------------------------
%u = fourierLS2(u,dx,(dt/2),umean,alpha,ep);            


plot_step = 1; 

index = 0; 
uold = u*0; 
inx_plot = 0; 

tol_level = 10; 
     
for j= 1:Nstep-1
    index = index + 1;
    time = time+dt; 
    
    % first linear half step ------------------------
    
    u = fourierLS2(u,dt);            

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
   
   
   
   u  = fourierLS2(u,dt);  
    
end

end


