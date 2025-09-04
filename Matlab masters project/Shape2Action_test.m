function [SA,SA1] = Shape2Action_test(Y,SD,tauD,sigmaD,alpha)
%function [SA,SARelArea] = Shape2Action(Y,SD,tauD,sigmaD,alpha)
% Energy difference between curved and flat state:
LD1=0.5*((Y(3,:)+sin(Y(2,:))./Y(1,:)-1).^2).*Y(1,:) + sigmaD*Y(1,:) + alpha*Y(3,:).*sin(Y(2,:)); % lagrangian dimensionless units
SA1=-trapz(SD,LD1) + tauD*Y(1,end); % energy dimensionless in curved state
%
LD=0.5*((Y(3,:)+sin(Y(2,:))./Y(1,:)-1).^2).*Y(1,:) + sigmaD*Y(1,:); % lagrangian dimensionless units
SA=-trapz(SD,LD) -alpha*(cos(Y(2,1))-cos(Y(2,end)))  + tauD*Y(1,end);
end