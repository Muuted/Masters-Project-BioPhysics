% Function finds the neck solution that matches an input alpha by varying psi2.
% Using DIMENSIONLESS UNITS
% Storing all solutions where kG/k=alpha_target as a vector of ode45 objects
% also storing all solutions with crossings t=tau.
%
%
function [Result]=ShapeAlpha(alpha_target,tauD,sigmaD)
try
    lambdaD=(0.5+sigmaD)^-0.5; % length scale asymptotic regime
    %r2D=2*tauD*lambdaD^2; % radius in state (2)
    % new procedure for finding radius r2 in initial state
    r2D=20; % using fixed r2D unless it is too small
    if r2D<0.8*tauD*lambdaD^2
        r2D=1.3*tauD*lambdaD^2;
    end
    %
    N1=4000; % number of psi2 values used for finding alpha in first step
    psi2_vector = -logspace(-0.5,-14,N1); % constructing the psi2 vector with logaritmic spacing.
    MaxTauCrossings=20; % maximum number of points with t=tau also used in Shape.
    MaxPsi2Zeros=20; % maximum number of psi2 values giving alpha=alpha_target
    alpha_matrix=zeros(MaxTauCrossings,N1); % initialize alpha vector
    h1 = waitbar(0,'Searching angles - first iteration');
    for m=1:N1
        waitbar(m/N1,h1)
        [~,alpha_matrix(:,m),~]=Shape(r2D,psi2_vector(m),tauD,sigmaD); % obtain alpha at a wide range of psi2 values
    end
    close(h1);
    alpha_deviation=(alpha_matrix-alpha_target)/alpha_target; % calculating fractional deviation in the kGs
    % finding psi2 values where alpha=alpha_target
    psi2_zeros=NaN(MaxTauCrossings,MaxPsi2Zeros);
    for k=1:MaxTauCrossings %iterate over solutions with t=tau
        Zeros=FindZeros1st(psi2_vector,alpha_deviation(k,:)); % find solutions with correct kG
        psi2_zeros(k,1:length(Zeros)) = Zeros; % putting psi2 values in matrix
    end
    %
    %% fine tuning psi2_zeros in second iteration:
    N2=3000; % no iterations in second step
    psi2_margin=0.1;
    maxlim=max(max(psi2_zeros))*(1-psi2_margin);
    minlim=min(min(psi2_zeros))*(1+psi2_margin);
    loglim2=log10(-maxlim);
    loglim1=log10(-minlim);
    psi2_vector = -logspace(loglim1,loglim2,N2); % new narrow psi2 vector
    % start second iteration
    alpha_matrix=zeros(MaxTauCrossings,N2); % initialize alpha vector
    h1 = waitbar(0,'Searching angles - second iteration');
    for n=1:N2
        waitbar(n/N2,h1)
        [~,alpha_matrix(:,n),~]=Shape(r2D,psi2_vector(n),tauD,sigmaD); % obtain alpha at a wide range of psi2 values
    end
    close(h1);
    alpha_deviation=(alpha_matrix-alpha_target)/alpha_target; % calculating fractional deviation in the kGs
    % finding psi2 values where alpha=alpha_target
    psi2_zeros=NaN(MaxTauCrossings,MaxPsi2Zeros);
    for k=1:MaxTauCrossings %iterate over solutions with t=tau
        Zeros=FindZeros2nd(psi2_vector,alpha_deviation(k,:)); % find solutions with correct kG
        psi2_zeros(k,1:length(Zeros)) = Zeros; % putting psi2 values in matrix
    end
    % plots for inspecting the alpha deviation curves as function of
    % psi2 and taus
    %
    % Code to plot result of shooting angles psi2 and detection of points
    % alpha_deviation=0
    %
    %
    % figure();
    % semilogx(psi2_vector,alpha_deviation(1,:),'.') %plotting deviation
    % hold on;
    % semilogx(psi2_vector,alpha_deviation(2,:),'.') %plotting deviation
    % semilogx(psi2_vector,alpha_deviation(3,:),'.') %plotting deviation
    % semilogx(psi2_vector,alpha_deviation(4,:),'.') %plotting deviation
    % semilogx(psi2_vector,alpha_deviation(5,:),'.') %plotting deviation
    % semilogx(psi2_vector,alpha_deviation(6,:),'.') %plotting deviation
    % line(xlim,[0 0],'Color','k');
    % plot(psi2_zeros,zeros(size(psi2_zeros)),'.','MarkerSize',12);
    % xlabel('start angle \psi_2 (rad)')
    % ylabel('Relative deviation \alpha')
    %%
    % We need to evaluate each of the solutions to give (Y,SD)
    NPoints=2000; % no of points for evaluating solution
    Result=struct([]);
    %
    for k=1:MaxTauCrossings %iterate over solutions with crossings t=tau
        for p=1:MaxPsi2Zeros % iterate psi2 that give alpha=alpha_target
            if ~isnan(psi2_zeros(k,p)) % pick the psi2zeros with non NaN values
                % Evaluate shape (k,p)
                [ShapeSolution,alpha_out,r0D]=Shape(r2D,psi2_zeros(k,p),tauD,sigmaD);
                SD=linspace(ShapeSolution.x(1),ShapeSolution.xe(k),NPoints); % make equally spaced S-points from s_start to s_end in the solution generated.
                Y=deval(ShapeSolution,SD);
                %
                if isempty(selfintersect(Y(1,:),Y(5,:))) % Remove Results with self intersections
                    % place all data relating to the combination TauCrossing (k)  and Psi2Zero (p) in structure array Result(k,p)
                    % Result structure is changing size for each iteration (not possible to preallocate)
                    Result(k,p).ShapeSolution=ShapeSolution;
                    Result(k,p).SD=SD;
                    Result(k,p).Y=Y;
                    Result(k,p).alpha_out=alpha_out(k); % extract only k value of alpha_out and r0D (corresponding to correct crossing):
                    Result(k,p).r0D=r0D(k);
                    % add hamiltonian
                    Result(k,p).HD=0.5*Y(1,:).*(Y(3,:).^2)-...
                        0.5*Y(1,:).*(sin(Y(2,:))./Y(1,:)-1).^2 +...
                        Y(4,:).*cos(Y(2,:))-sigmaD*Y(1,:);
                    % place input variables in Result
                    Result(k,p).r2D=r2D; % radius in start point
                    Result(k,p).psi2=psi2_zeros(k,p); % psi2
                    Result(k,p).tauD=tauD; % line tension
                    Result(k,p).sigmaD=sigmaD; % tension
                    Result(k,p).psi2_vector=psi2_vector; % psi2 vector used for shooting
                    Result(k,p).alpha_deviation=alpha_deviation; % relative deviation to alpha target
                    Result(k,p).psi2_zeros=psi2_zeros; % psi2
                    [Result(k,p).SA,Result(k,p).SA_flat,Result(k,p).DeltaSA]=Shape2Action(Result(k,p).Y,Result(k,p).SD,tauD,sigmaD,alpha_target);
                    Result(k,p).Area=-Y(6,end);
                    Result(k,p).ExcessArea=Result(k,p).Area-pi*r2D.^2;
                else
                    disp('self intersection detected and removed');
                    continue
                end
            else
                continue
            end
        end
    end
catch
    disp('An error was generated');
    Result=struct([]);
    if  isgraphics(h1)
        close(h1)
    end
end
end
