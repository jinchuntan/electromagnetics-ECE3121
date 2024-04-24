function Efieldplot(opt)
    % Written by Tan Jin Chun
    % Last Modified: 4/9/2023
    % Generation of Electric Field
    % Function Efieldplot produces a variety
    % of plots for collections of line charges.
    % Charges and their locations are set up below.

    % Modelling the 4 electrodes (Mode = Central Pour)
    Ql(1)=10; x(1)=0; y(1)=1.0; % Electrode A
    Ql(2)=10; x(2)=0; y(2)=-1.0; % Electrode C
    Ql(3)=10; x(3)=1.0; y(3)=0; % Electrode B
    Ql(4)=10; x(4)=-1.0; y(4)=0; % Electrode D
    Ql(5)=0; x(5)=0; y(5)=0;
    Ql(6)=0; x(6)=0; y(6)=0;

    % Set up grid.
    [X,Y] = meshgrid(-2:0.1:2);

    % Set maximum number of line charges
    nMAX=6;

    % Pre-assign a number of arrays
    Ex=zeros(length(X),length(Y),nMAX);
    Ey=zeros(length(X),length(Y),nMAX);
    EX=zeros(length(X),length(Y));
    EY=zeros(length(X),length(Y));
    rhosq=ones(length(X),length(Y),nMAX);

    for n = 1:length(Ql)
        if Ql(n)==0
            nmax=n-1; break;
        elseif n >= nMAX
            % No of charges cannot exceed nMAX
            nmax=nMAX;
        end
    end

    for n = 1:nmax
        rhosq(:,:,n) = (X-x(n)).^2 + (Y-y(n)).^2;
        Ex(:,:,n) = Ql(n)*(X-x(n))./(rhosq(:,:,n)+eps);
        Ey(:,:,n) = Ql(n)*(Y-y(n))./(rhosq(:,:,n)+eps);

        for i=1:length(X)
            for j=1:length(Y)
                if rhosq(i,j,n) < 0.005
                    Ex(i,j,n) = 0;
                    Ey(i,j,n) = 0;
                end
            end
        end

        EX = EX + Ex(:,:,n);
        EY = EY + Ey(:,:,n);
    end

    EX = EX/(2*pi);
    EY = EY/(2*pi);

    % Generate plots
    figure(1)
    quiver(X,Y,EX,EY)
    axis square
    xlabel('X'); ylabel('Y');
    title('Electric Field Vectors');

    figure(2)
    plot(X(1+(length(X)-1)/2,:),EX(:,1+(length(X)-1)/2), ...
        X(1+(length(X)-1)/2,:),EY(:,1+(length(X)-1)/2))
    xlabel('Y'); ylabel('E field');
    legend('EX','EY')
    title('Electric Field Components on Y axis');

    figure(3)
    plot(Y(:,1+(length(X)-1)/2),EX(1+(length(X)-1)/2,:), ...
        Y(:,1+(length(X)-1)/2),EY(1+(length(X)-1)/2,:))
    xlabel('X'); ylabel('E field');
    legend('EX','EY')
    title('Electric Field Components on X axis');

    if nargin == 1
        if strcmp(opt,'nosections')
            return;
        end
    end

    figure(4)
    plot(X(1,:),EX')
    legend(num2str(Y(:,1)), 'location', 'BO')
    grid on
    title('Electric Field Component for various Y')
    xlabel('X'); ylabel('X component of E field');

    figure(5)
    plot(X(1,:),EY')
    legend(num2str(Y(:,1)), 'location', 'BO')
    grid on
    title('Electric Field Component for various Y')
    xlabel('X'); ylabel('Y component of E field');
end


