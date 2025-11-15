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
ymin = 2.0;
ymax = 2.025;
xmin = 39;
xmax = 40;
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
clc
figure(1)
plot_color=jet(N_tauD);

for n=1:N_tauD
    plot(ExcessAreaCurved{n},NeckRadiusCurved{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
    hold on
    %flat_area_point(n) = 
    plot(ExcessAreaFlat{n},NeckRadiusFlat{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
end
for n=1:N_tauD
    N = size(ExcessAreaCurved{n});
    N(2)
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
            if ymin < r0D && r0D < ymax
                if xmin < AE && AE < xmax
                    AE;
                    r0D;
                    plot(AE,r0D,"*","MarkerSize",15)
                    tau = tau_ref;
                    sigma = sigma_ref;
                    psi_L = psiL_ref;
                end
            end
        end
    end
end
xlim([-60 50])
ylim([0 7])

%%
tau
sigma
psi_L

%% Saving the data to .txt file for Python
clc
name_list = [
    "sigma curved"
    "tau curved"
    "psi2 curved" 
    "A excess curved"
    "dSA curved"
    "Neck r0 curved" 
    "r1 curved" 
    "sigma flat"
    "tau flat"
    "psi2 flat"
    "A excess flat"
    "dSA flat"
    "Neck r0 flat"
    "r1 flat"
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
data_path = "C:\Users\AdamSkovbjergKnudsen\Documents\GitHub\Masters-Project-BioPhysics\Surface model\2D sim results\";
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
save_cell = cell(2,2);
save_cell{2,1} = "ExcessAreaCurved";
save_cell{2,2} = save_matrix_curved(4,:);

save_cell{3,1} = "NeckRadiusCurved";
save_cell{3,2} = save_matrix_curved(6,1:100);

save_cell
writecell(save_cell,strcat(data_path,"ExcessAreaCurved.dat"));
writecell(save_cell,strcat(data_path,"ExcessAreaCurved.txt"));

%%

save_cell = cell(2,2);
save_cell{1,1} = "index";
len = length(save_matrix_curved(6,:));
save_cell{1,2} = linspace(1,len,len);
save_cell{2,1} = "NeckRadiusCurved";
save_cell{2,2} = save_matrix_curved(6,:);
save_cell
writecell(save_cell,strcat(data_path,"NeckRadiusCurved.dat"));



%%
writecell(ExcessAreaCurved,strcat(data_path,"ExcessAreaCurved.dat"));
writecell(NeckRadiusCurved,strcat(data_path,"NeckRadiusCurved.dat"));

writecell(ExcessAreaCurved,strcat(data_path,"ExcessAreaCurved.txt"));
writecell(NeckRadiusCurved,strcat(data_path,"NeckRadiusCurved.txt"));


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