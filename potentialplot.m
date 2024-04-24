function potentialplot
% Written by Tan Jin Chun
% Last Modified: 4/9/2023

% Set up grid.
[X,Y] = meshgrid(-2:.01:2);

% Set maximum number of line charges
nMAX=6;

% Pre-assign a number of arrays
x=zeros(1,nMAX);y=zeros(1,nMAX);Q=zeros(1,nMAX);

% Read input data
inputfile

% Entering a zero charge or omitting a charge will result in charges with higher index numbers
% being neglected.
for n =1:length(Ql)
    if Ql(n)==0
        nmax=n-1;break
    elseif n >= nMAX

        % No of charges cannot exceed nMAX
        nmax=nMAX;
    end
end

% The preceeding code does not cope with all cases of missing charges satisfactorily.
% Pre-assign more arrays
v=zeros(length(X),length(Y),nmax);
V=zeros(length(X),length(Y));
rhosq=ones(length(X),length(Y),nMAX);

% For charge number n, calculate field
for n =1:nmax
    rhosq(:,:,n) = (X-x(n)).^2 + (Y-y(n)).^2;
v(:,:,n) = - Ql(n)*log(rhosq(:,:,n)+eps)/2;

% To avoid trouble with infinite values
% set field to zero near charges
for i=1:length(X)
    for j=1:length(Y)
        if rhosq(i,j,n) < 0.005
        v(i,j,n) = 0;
        end
    end
end

% Sum contributions to field
V = V + v(:,:,n);

end

% Scale fields. Note epsilon is omitted.
V = V/(2*pi);

% plot showing equipotentials
figure(1)
[C h]=contour(X,Y,V);
clabel(C,h);
axis square
grid on
