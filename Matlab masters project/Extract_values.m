clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}) % set to real black
%set(groot,'defaultAxesFontName','Helvetica')
%set(groot,'DefaultLegendFontName','Helvetica')
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%% Load the data

file = "MS_tauD_1_to_7_NtauD_20_sigmaD_m0p4_to_2_NsigmaD_80.mat";

data_set = load(file);
clc
MS = data_set.MS;
[N_tauD,N_sigmaD]=size(MS);

psi1_cut=-pi;
sigma_list_Curved = cell(1,N_tauD);
sigma_list_Flat = cell(1,N_tauD);

tau_list_Curved = cell(1,N_tauD);
tau_list_Flat = cell(1,N_tauD);

psi_L_list_Curved = cell(1,N_tauD);
psi_L_list_Flat = cell(1,N_tauD);

ExcessAreaCurved=cell(1,N_tauD);
ExcessAreaFlat=cell(1,N_tauD);

DeltaSACurved=cell(1,N_tauD);
DeltaSAFlat=cell(1,N_tauD);

NeckRadiusCurved=cell(1,N_tauD);
NeckRadiusFlat=cell(1,N_tauD);

r1Curved=cell(1,N_tauD);
r1Flat=cell(1,N_tauD);



%%
for n=1:N_tauD
    ah=1; % index to high curvature points
    al=1; % index to low curvature points
    for m=1:N_sigmaD % loop sigmaD
        [R,C]=size(MS(n,m).Result);
        if R~=0 && C~=0 % proceed if Result is not empty
            for k=1:R % loop over crossings t=tau
                for p=1:C % loop over psi2zero values
                    if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty                       
                        if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut && MS(n,m).Result(k,p).ExcessArea>0 % if high curved and excessarea>0
                           % extracting data points for high curved necks
                            ExcessAreaCurved{n}(ah)=MS(n,m).Result(k,p).ExcessArea; % Excess area of curved neck
                            DeltaSACurved{n}(ah)=MS(n,m).Result(k,p).DeltaSA; % relative energy of curved neck
                            NeckRadiusCurved{n}(ah)=min(MS(n,m).Result(k,p).Y(1,:)); % neck radius as minimum radius of shape
                            r1Curved{n}(ah)=MS(n,m).Result(k,p).Y(1,end); % radius at free edge

                            sigma_list_Curved{n}(ah) = MS(n,m).Result(k,p).sigmaD;
                            tau_list_Curved{n}(ah) = MS(n,m).Result(k,p).tauD;
                            psi_L_list_Curved{n}(ah) = MS(n,m).Result(k,p).psi2;
                            ah=ah+1;
                        elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                            % extracting data points for low curved necks
                            ExcessAreaFlat{n}(al)=MS(n,m).Result(k,p).ExcessArea;
                            DeltaSAFlat{n}(al)=MS(n,m).Result(k,p).DeltaSA;
                            NeckRadiusFlat{n}(al)=min(MS(n,m).Result(k,p).Y(1,:)); % neck radius as minimum radius of shape
                            r1Flat{n}(al)=MS(n,m).Result(k,p).Y(1,end); % radius at free edge

                            sigma_list_Flat{n}(al) = MS(n,m).Result(k,p).sigmaD;
                            tau_list_Flat{n}(al) = MS(n,m).Result(k,p).tauD;
                            psi_L_list_Flat{n}(al) = MS(n,m).Result(k,p).psi2;

                            al=al+1;
                            legend_text{n}=['$\tilde{\tau}$ = ' num2str(MS(n,m).Result(k,p).tauD,2)]; % construct legend text
                        end
                    end
                end
            end
        end
    end
    % removing degenerate areapoints with higher energy (flat necks)
    id_ok_flat=find(diff(ExcessAreaFlat{n})>0);
    ExcessAreaFlat{n}=ExcessAreaFlat{n}(id_ok_flat);
    DeltaSAFlat{n}=DeltaSAFlat{n}(id_ok_flat);
    NeckRadiusFlat{n}=NeckRadiusFlat{n}(id_ok_flat);
    r1Flat{n}=r1Flat{n}(id_ok_flat);
end
%%
%% span of the point

%%
close all
clc
figure(1)
hold on
grid on
xlim([0 40])
ylim([0 2.25])
plot_color=jet(N_tauD);
xmin_list = [3 4 4.8 5.2 5.8 6.5 7.4 9.4 10.6 12 13.4 15 16.6 18 20 21.5 23.5 25 27 29 31 33.5 36.2 39.2];%xmin
xmax_list = [4 4.1 4.82 5.5 6 6.7 7.42 9.5 10.7 12.1 13.5 15.1 16.7 18.5 20.5 22 24 26 27.5 29.75 32 34 36.4 39.6];%xmax
ymin_list = [0.1 0.145 0.194  0.2 0.24 0.3 0.36 0.47 0.56 0.64 0.73 0.82 0.91 1 1.1 1.2 1.3 1.4 1.45 1.55 1.65 1.75 1.9 2.02];%ymin
ymax_list = [0.15 0.148 0.2 0.24 0.28 0.32 0.364 0.5 0.57 0.65 0.74 0.83 0.93 1.1 1.2 1.3 1.35 1.45 1.525 1.6 1.725 1.8 1.91 2.03];%ymax


tau_return_list = [];
sigma_return_list = [];
psi2_return_list = [];

for n=1:N_tauD
    plot(ExcessAreaCurved{n},NeckRadiusCurved{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
    hold on
    %flat_area_point(n) = 
    plot(ExcessAreaFlat{n},NeckRadiusFlat{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
end

pause(2)

for i=1:length(ymin_list)
    for n=1:N_tauD
        N = size(ExcessAreaCurved{n});
        for k=1:N(2)
            for m=1:2
                if m==1
                    AE = ExcessAreaCurved{n}(k);
                    r0D= NeckRadiusCurved{n}(k);

                    sigma_ref = sigma_list_Curved{n}(k);
                    tau_ref = tau_list_Curved{n}(k);
                    psiL_ref = psi_L_list_Curved{n}(k);
                end
                if m==2
                    AE = ExcessAreaFlat{n}(k);
                    r0D= NeckRadiusFlat{n}(k);

                    sigma_ref = sigma_list_Flat{n}(k);
                    tau_ref = tau_list_Flat{n}(k);
                    psiL_ref = psi_L_list_Flat{n}(k);
                end
                if ymin_list(i) < r0D && r0D < ymax_list(i)
                    if xmin_list(i) < AE && AE < xmax_list(i)
                        AE;
                        r0D;
                        plot(AE,r0D,"*","MarkerSize",15)
                        tau = tau_ref;
                        sigma = sigma_ref;
                        psi_L = psiL_ref;
                        tau_return_list(i) = tau_ref;
                        sigma_return_list(i) = sigma_ref;
                        psi2_return_list(i) = psiL_ref;
                    end
                end
            end
        end
    end
    pause(0.5)
end


%%
clc
[tau_return_list]
%%
clc
[sigma_return_list]
%%
clc
[psi2_return_list]

%% Saving the data to .txt file for Python
clc
name_list =[
    "sigma_list_Curved"
    "tau_list_Curved"
    "psi_L_list_Curved"
    "ExcessAreaCurved"
    "DeltaSACurved"
    "NeckRadiusCurved"
    "r1Curved"
    "sigma_list_Flat"
    "tau_list_Flat"
    "psi_L_list_Flat"
    "ExcessAreaFlat"
    "DeltaSAFlat"
    "NeckRadiusFlat"
    "r1Flat"
    ];
data_list =[
    sigma_list_Curved
    tau_list_Curved
    psi_L_list_Curved
    ExcessAreaCurved
    DeltaSACurved
    NeckRadiusCurved
    r1Curved
    sigma_list_Flat
    tau_list_Flat
    psi_L_list_Flat
    ExcessAreaFlat
    DeltaSAFlat
    NeckRadiusFlat
    r1Flat
];

%%
data_path = "C:\Users\AdamSkovbjergKnudsen\Documents\GitHub\Masters-Project-BioPhysics\Surface model\2D sim results\matlab data transfer files\";
save_name = name_list(1);
save_data_for_python(data_path,save_name,sigma_list_Curved,".dat")

save_name = name_list(2);
save_data_for_python(data_path,save_name,tau_list_Curved,".dat")


save_name = name_list(3);
save_data_for_python(data_path,save_name,psi_L_list_Curved,".dat")


save_name = name_list(4);
save_data_for_python(data_path,save_name,    ExcessAreaCurved,".dat")

save_name = name_list(5);
save_data_for_python(data_path,save_name,    DeltaSACurved,".dat")

save_name = name_list(6);
save_data_for_python(data_path,save_name,    NeckRadiusCurved,".dat")

save_name = name_list(7);
save_data_for_python(data_path,save_name,    r1Curved,".dat")




save_name = name_list(8);
save_data_for_python(data_path,save_name,sigma_list_Flat,".dat")

save_name = name_list(9);
save_data_for_python(data_path,save_name,tau_list_Flat,".dat")


save_name = name_list(10);
save_data_for_python(data_path,save_name,psi_L_list_Flat,".dat")


save_name = name_list(11);
save_data_for_python(data_path,save_name,    ExcessAreaFlat,".dat")

save_name = name_list(12);
save_data_for_python(data_path,save_name,    DeltaSAFlat,".dat")

save_name = name_list(13);
save_data_for_python(data_path,save_name,    NeckRadiusFlat,".dat")

save_name = name_list(14);
save_data_for_python(data_path,save_name,    r1Flat,".dat")



%%


function [] = save_data_for_python(data_path,save_name,data,file_type)
    N_tauD = length(data);
    save_cell = cell(N_tauD,2);
    max_len = 0;
    for n=1:N_tauD
        M = [];
        len = length(data{n});
        if max_len < len
            max_len = len;
        end
        for j=1:len
            M(j) = data{n}(j);
        end
        save_cell{n+1,1} = strcat(save_name," ","n=",string(n));
        save_cell{n+1,2} = M;
    end
    save_cell{1,1} = "index";
    save_cell{1,2} = linspace(1,max_len,max_len);
    writecell(save_cell,strcat(data_path,save_name,file_type));
end