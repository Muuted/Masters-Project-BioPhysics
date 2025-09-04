function [SA,SA_flat,DeltaSA] = Shape2Action(Y,SD,tauD,sigmaD,alpha)
%function [SA,SARelArea] = Shape2Action(Y,SD,tauD,sigmaD,alpha)
% Energy difference between curved and flat state:
%
LD=0.5*((Y(3,:)+sin(Y(2,:))./Y(1,:)-1).^2).*Y(1,:) ; % lagrangian dimensionless units - without sigmaD (17/05/23). No difference if solutions with equal areas are compared.
SA=-trapz(SD,LD) -alpha*(cos(Y(2,1))-cos(Y(2,end)))  + tauD*Y(1,end); % Action. Gaussian term is boundary term.
SA_flat=-0.5*(Y(6,end)/(2*pi)); % flat reference state without a hole
DeltaSA=SA-SA_flat; % Energy difference to a flat reference state without a hole
%
end