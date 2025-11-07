clc;
clearvars;
clear all;


%% My simulation units variables
c0 = 0.25;%1.0;%.25;
k  = 1.;
sigma = k*(c0^2); 
tau = 1;%0.1;

%% The  dimless units conversions.
lc = 1/c0;
sigma_c = k*(c0^2);
tau_c = k*c0;
c0_c = c0;
k_c = k;

%% The dimless units
clc
"the dimless units"
k = k/k_c;
tauD= 1%tau/tau_c %1
sigmaD= 0.1%sigma/sigma_c ;%0.1
c0 = c0/c0_c;
r2D = 20;%tau*1.01/(k*c0^2)%20
psi2 = -7.3648e-8; % -2.8362e-8

%% Finding the solution

[ShapeSolution,alpha,r0D]=Shape(r2D,psi2,tauD,sigmaD);

m = 0;
for i= 1:length(ShapeSolution.y(1,:))
    if ShapeSolution.y(4,i) < tauD
        m = i;
        break
    end
end
m
r = ShapeSolution.y(1,1:m);
z = ShapeSolution.y(5,1:m);
psi = ShapeSolution.y(2,1:m);
%r = ShapeSolution.y(1,:);
%z = ShapeSolution.y(5,:);

figure(1)
plot(r,z)
%%
figure(2)
plot(psi)
%plot(Result(k,p).Y(1,:),Result(k,p).Y(5,:)) % plot (r,z) in solution
%% Making new matrix
m = length(ShapeSolution.y(1,:))
save_matrix = zeros(6,m);
size(save_matrix)
for i=2:6
    for j=1:m
        save_matrix(i,j) = ShapeSolution.y(i-1,j);
    end
end
%% Saving integration results

data_path = "C:\Users\AdamSkovbjergKnudsen\Documents\GitHub\Masters-Project-BioPhysics\Matlab masters project\saved data\";
name = "Compare integration results";
save_file_name = strcat(data_path,name,".txt");
writematrix(save_matrix,save_file_name );
%writematrix(ShapeSolution.y,save_file_name );
"done"
