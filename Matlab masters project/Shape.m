% This function outputs alpha and r0D when finding the solution. To be used
% in a routine that adjusts psi2 to obtain alpha.
% Using DIMENSIONLESS UNITS
% Revision Okt 2022: Not stopping integration when t=tau, but detecting all
% crossings!
function [ShapeSolution,alpha,r0D]=Shape(r2D,psi2,tauD,sigmaD)
%% calculated quantities
sigmaD;
lambdaD=(0.5+sigmaD)^(-0.5); % length scale in assymptotic solution
p2D=r2D/lambdaD; % variable in bessel function
AD=psi2/besselk(1,p2D); % modified bessel function of second kind, order=1.
psidot2D=(AD/lambdaD)*(-besselk(0,p2D)-besselk(1,p2D)./p2D); % Derivative of psi2
t2D=(sigmaD+0.5)*r2D; % t2 found from H=0 when psi=0 and psidot=0
z2D=0; % offset z axis
Area0=0; % start area for integration
%% Solve differential equations:
% Define the range of s values we want to solve for (unit of s is meters):
srangeD=[0 -r2D-10]; % (meters). Note that we integrate in negative s direction from s2=0 to s1<0.
%srangeD = linspace(0, -r2D-10,5000);
% Integration from s2 to s1 is necessary because the state s2 is known (asymptotic), but the state s1 is not.
initialvalues=[... % initial values written into a vector
    r2D ...    % r(s2)
    psi2 ...  % psi(s2)
    psidot2D ...% psidot(s2)
    t2D ...   % t(s2)
    z2D...   % z(s2)
    Area0...   % Area(s2)
    ];
% Next we call the MATLAB function ode45 that solves differential equations numerically.
% The equations must be a set of coupled (ordinary) first order differential
% equations. Note that ode45 uses the syntax @functionname which is a
% function handle (reference to a function).
%
%Opt    = odeset('Events', @EdgeTensionReached); % event criterium for integration is that t=tau
Opt    = odeset('Events', @EdgeTensionReached,'RelTol', 1E-6); % event criterium for integration is that t=tau
%Opt    = odeset('Events', @EdgeTensionReached,'RelTol', 1E-12,'AbsTol',1E-12); % event criterium for integration is that t=tau
%Opt    = odeset('Events', @EdgeTensionReached); % event criterium for integration is that t=tau
%
%Opt    = odeset('Events', @EdgeTensionReached); % event criterium for integration is that t=tau
% Event function that detects when t=tau during integration
    function [value, isterminal, direction] = EdgeTensionReached(~, Y)
        value      = Y(4) - tauD;
        isterminal = 0;   % 0=Do not stop the integration when Y(4)=tau
        direction  = 0;
    end
%
ShapeSolution=ode45(@NeckDiffEquations,srangeD,initialvalues,Opt); % obtain solution of shape equations. Includes all solutions t=tau in one object.
% Evaluate kG/alpha and r0D at all crossings t=tau
MaxTauCrossings=20; % maximum number of times that t=tau during integration
alpha=NaN(MaxTauCrossings,1); % initialize vector holding alpha. Needs to be fixed size to plot for varying psi2
r0D=NaN(MaxTauCrossings,1); % initialize r0D vector
for k=1:length(ShapeSolution.ie) % iterate over crossings t=tau
    Y_end=deval(ShapeSolution,ShapeSolution.xe(k)); % evaluate solution vector at crossing #k
    r1D=Y_end(1); % radius at crossing #k
    psi1=Y_end(2); % angle at crossing #k
    psidot1D=Y_end(3); % dpsi/ds at crossing #k
    alpha(k)=(1-psidot1D)*r1D/sin(psi1)-1; %alpha=kG/k at crossing #k
    % Finding initial hole radius r0D based on integrated membrane area
    Area=abs(Y_end(6)); % area of curved state. Use absolute value because result is negative.
    r0D(k)=sqrt((pi*r2D^2 - Area)/pi); % find r0D by setting Area= Area_flat
end
%
%% Define differential equations:
    function dy = NeckDiffEquations(~,y)
        % This function defines the set of differential equations dy defining the membrane shape around a cirular hole.
        % Equations must follow the format: dy = f(s,y), where dy and y can be vectors and dy=dy/ds
        % Here s is the independent variable. The dependent variables y(1),y(2).. will be found when solving using ode45.
        % We define the set of differential equations as a vector = dy, where:
        % y(1)=rD;    dy(1)=drD/dsD
        % y(2)=psi;  dy(2)=dpsi/dsD
        % y(3)=dpsi/dsD; dy(3)=d2psi/dsD2
        % y(4)=tD;    dy(4)=dtD/dsD
        % y(5)=zD;    dy(5)=dzD/dsD
        % y(6)=Area;  dy(6)=dArea/ds
        dy=[...
            cos(y(2));... % dy(1)=dr/ds
            y(3);... % dy(2)=dpsi/ds
            -y(3).*cos(y(2))./y(1) + sin(y(2)).*cos(y(2))./y(1).^2 + y(4).*sin(y(2))./(y(1));...% dy(3)=d2psi/ds2
            0.5*( (y(3)-1).^2 - (sin(y(2))./y(1)).^2 ) + sigmaD;... % dy(4)=dt/ds
            sin(y(2));... % dy(5)=dz/ds
            2*pi*y(1) % dy(6)=dArea/ds
            ];
    end
end