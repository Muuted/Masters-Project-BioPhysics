% This script is for testing the neck shape and plotting the result
% Using DIMENSIONLESS UNITS
%% input parameters
clear all;
%
alpha_target=-0.75; % alpha=kG/k
tauD=1; % linetension
sigmaD=-0.2; % tension
lambdaD=(0.5+sigmaD)^-0.5; % length scale asymptotic regime
r2D=2*tauD*lambdaD^2; % radius in (2)
psi2=-0.0030;
%
%% 
[ShapeSolution,alpha,r0D]=Shape(r2D,psi2,tauD,sigmaD);
%
NPoints=1000; % number of points when evaluating solution with deval
%SD=linspace(ShapeSolution.x(1),ShapeSolution.xe,NPoints); % make equally spaced S-points from s_start to s_end in the solution generated.
SD=linspace(ShapeSolution.x(1),ShapeSolution.x(end),NPoints); % make equally spaced S-points over whole S-range.
Y=deval(ShapeSolution,SD); % evaluate solution from structure using S-vector above.
% Y(1)=rD;  
% Y(2)=psi;  
% Y(3)=psidotD; 
% Y(4)=tD;  
% Y(5)=zD;    
% Hamiltonian :
HD=0.5*Y(1,:).*(Y(3,:).^2) ...
    -0.5*Y(1,:).*(sin(Y(2,:))./Y(1,:)-1).^2 ...
    +Y(4,:).*cos(Y(2,:)) ...
    -sigmaD*Y(1,:);
% Lagrangian, including gaussian curvature in Lagrangian for simplicity in
% integration:
% LD=0.5*((Y(3,:)+sin(Y(2,:))./Y(1,:)-1).^2).*Y(1,:) + sigmaD*Y(1,:) + alpha*Y(3,:).*sin(Y(2,:)); % lagrangian dimensionless units
% GpD=-trapz(SD,LD) + tauD*Y(1,end); % energy dimensionless in curved state
% rFlatD=r2D+SD; % radius vector in flat state
% L0D=0.5*rFlatD + sigmaD*rFlatD; %
% Gp0D=-trapz(SD,L0D) + tauD*r0D; % action in flat state
% DeltaGD=GpD-Gp0D; % energy difference 
%
%% plotting
figure();
%**** plot profile (r,z)
subplot(2,2,1)
plot(Y(1,:),Y(5,:)) % plot (r,z) in solution
hold on;
plot(ShapeSolution.ye(1,:),ShapeSolution.ye(5,:),'og');
title('Neck profile')
xlabel('rD');
ylabel('zD');
axis equal
%***********plot tension (s,t) *************
subplot(2,2,2)
plot(SD,Y(4,:),'.k-') % plot t as function of s
hold on
axis equal
plot(ShapeSolution.xe,ShapeSolution.ye(4),'og')
xlabel('sD ');
ylabel('tD');
title('Tension versus S')
%******** plot t as function of r
subplot(2,2,3)
plot(Y(1,:),Y(4,:),'.k-') % plot t as function of r
hold on
axis equal
tas=(0.5+sigmaD)*Y(1,:);
plot(Y(1,:),tas,'--k');
plot(ShapeSolution.ye(1,:),ShapeSolution.ye(4,:),'og')
xlabel('rD ');
ylabel('tD');
title('Tension versus S')
%

% %************plot (r,z) symmetric around hole ***********
% subplot(2,2,3)
% plot(Y(1,:),Y(5,:)) % profile
% hold on
% axis equal
% plot(-Y(1,:),Y(5,:)) % plot profile other side 
% plot(r0D,0,'.r');
% plot(-r0D,0,'.r');
% xlabel('rD ');
% ylabel('zD');
% title('Neck shape')
% %************* plot curvatures
% C1=Y(3,:);
% C2=sin(Y(2,:))./Y(1,:);
% SumC=C1+C2;
% subplot(2,2,4)
% plot(SD,C1,'r'); % profile
% hold on
% plot(SD,C2,'b');
% plot(SD,SumC,'g')
% xlabel('sD');
% ylabel('Curvature (m^{-1})');
% title('Curvatures')


