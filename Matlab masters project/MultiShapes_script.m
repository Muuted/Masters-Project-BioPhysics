%% This script is for generating multiple solutions to the shape equations by varying tauD and sigmaD in a matrix
% Output of script is the structure MS containing one Result for each
% (tauD,sigmaD) combination
% 
clearvars;
warning('off','all')
% input parameters
alpha_target=-0.75;
% Specify line tension values taudD (input):
tauD_min=0.2;
tauD_max=2;
N_tauD=20; % number of line tension values between tauD_min and tauD_max
tauD=linspace(tauD_min,tauD_max,N_tauD);
% Specify tension values (input):
sigmaD_min=-0.4;
sigmaD_max=2;
N_sigmaD=80;
sigmaD=linspace(sigmaD_min,sigmaD_max,N_sigmaD);
%
MS=struct([]); % initialize output structure
for n=1:N_tauD% loop tauD
    for m=1:N_sigmaD % loop sigmaD
        disp(n),disp(m);
        MS(n,m).Result=ShapeAlpha(alpha_target,tauD(n),sigmaD(m));
    end
end
