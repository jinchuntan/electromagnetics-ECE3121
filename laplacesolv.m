function [V,it,error]=laplacesolv(maxit,maxerror)
% Written by Tan Jin Chun
% Last Modified: 4/9/2023 

% Function [V,it,error]=laplacesolv
% This function solves for the potential distribution within the square domain using the successive over-relaxation method.
% The potential is iteratively updated until either a maximum number of iterations is reached or the error in the solution is less than a specified value.
% At the end of the function, the computed potential distribution V is visualized using contour plots and a 3D surface plot.

% Solves Laplace's equation for a 2 D square region by the difference equation method.

tic
% Some variables are declared global so that
% they will be available from the function setVBs
global Va Vb Vc Vd;
global X Y Vinitdisp;
Vinit = flipud(Vinitdisp);
V = Vinit;

% set default values
% Changing the number of iterations 
itmaxdef = 5;
errormaxdef = 1e-10;

if nargin==1
    itmax=maxit;
else
    itmax=itmaxdef;
end

if nargin==2
    errormax=maxerror;
    itmax=maxit;
else
    errormax=errormaxdef;
end

for it = 1:itmax
    Vold=V;
    for i=2:length(X)-1
        for j=2:length(Y)-1
            if  Vinit(i,j)==0
                V(i,j)=(V(i-1,j)+V(i,j-1)+V(i+1,j)+V(i,j+1))/4;
            else
                V(i,j)=V(i,j);
            end
        end
    end
    error=max( max(abs(Vold-V)) )/max( max(abs(Vinit)) );
    if error < errormax, break, end
end
it     %for diagnostic use
error  %for diagnostic use
Vold-V;


figure(6) % plot showing equipotentials
[C h]=contour(X,Y,V);
clabel(C,h);
text(-0.2,-1.9,['Va = ',num2str(Va)]);
text(1.7,0,['Vb = ',num2str(Vb)]);
text(-0.2,1.9,['Vc = ',num2str(Vc)]);
text(-1.9,0,['Vd = ',num2str(Vd)]);
axis square
grid on
xlabel('X');ylabel('Y')
title('The Equipotentials of the 4 Electrodes')
V;
toc

%% Adding few lines of code here
figure(7)
surf(X,Y,V)
title('The plot of the equipotential surface of the 4 electrodes')