%% Script to plot the energy difference DeltaSA for the neck versus the neck area.
% For figure 5 + 6 in revised manuscript march 2024
%
% Setting the default properties of the axis
clear all
close all
clc

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}) % set to real black
%set(groot,'defaultAxesFontName','Helvetica')
%set(groot,'DefaultLegendFontName','Helvetica')
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%% Load the data
data_set = load("MS_tauD_1_to_7_NtauD_20_sigmaD_m0p4_to_2_NsigmaD_80.mat");
%%

clc
MS = data_set.MS
size(MS)
n = 4;
m = 2;
Results = MS(n,m).Result(m)

Results.psi2

%%

Results()
%% plot (Area,DeltaSA) for fixed tauD. 
% FIGURE 1

[N_tauD,N_sigmaD]=size(MS);
psi1_cut=-pi; % cutoff between large and small type necks
f1=figure();
subplot(1,3,2);
plot_color=jet(N_tauD);
%plot_color='b';
ExcessAreaCurved=cell(1,N_tauD);
DeltaSACurved=cell(1,N_tauD);
ExcessAreaFlat=cell(1,N_tauD);
DeltaSAFlat=cell(1,N_tauD);
NeckRadiusCurved=cell(1,N_tauD);
NeckRadiusFlat=cell(1,N_tauD);
r1Curved=cell(1,N_tauD);
r1Flat=cell(1,N_tauD);
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
                            ah=ah+1;
                        elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                            % extracting data points for low curved necks
                            ExcessAreaFlat{n}(al)=MS(n,m).Result(k,p).ExcessArea;
                            DeltaSAFlat{n}(al)=MS(n,m).Result(k,p).DeltaSA;
                            NeckRadiusFlat{n}(al)=min(MS(n,m).Result(k,p).Y(1,:)); % neck radius as minimum radius of shape
                            r1Flat{n}(al)=MS(n,m).Result(k,p).Y(1,end); % radius at free edge
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
% plotting extracted data points
for n=1:N_tauD
    plot(ExcessAreaCurved{n},DeltaSACurved{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
    hold on
    flatSApoint(n)=plot(ExcessAreaFlat{n},DeltaSAFlat{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
end
h=gca;
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5); % Vertical line at zero excess area
legend(flatSApoint,legend_text); % add new legend for each value of tauD
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$','Interpreter','latex')
ylabel('Relative energy  $\Delta \tilde{S}=\tilde{S}_{neck}-\tilde{S}_{flat}$','Interpreter','latex')
title('Neck energy versus area')
%% plotting selected neck profiles for figure in paper
n=2;
delta=-3;
offset=0;
num_sigma_flat=1;
num_sigma_curved=1;
sigma_values=[6 9 12 15 18 30 50];
sigma_text_flat=cell(1,1);
sigma_text_curved=cell(1,1);
for m=sigma_values % loop sigmaD    
    [R,C]=size(MS(n,m).Result);
    if R~=0 && C~=0 % proceed if Result is not empty
        for k=1:R % loop over crossings t=tau
            for p=1:C % loop over psi2zero values
                if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty
                    if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut && MS(n,m).Result(k,p).ExcessArea>0 % if high curved and excessarea>0
                        % extracting data points for high curved necks
                        subplot(1,3,2)
                        plot(MS(n,m).Result(k,p).ExcessArea,MS(n,m).Result(k,p).DeltaSA,'.k','MarkerSize',13,'HandleVisibility','off') % mark point
                        hold on
                        sigma_text_curved{num_sigma_curved}=['$\tilde{\sigma}$ = ' num2str(MS(n,m).Result(k,p).sigmaD,2)];
                        num_sigma_curved=num_sigma_curved+1;                        
                        subplot(1,3,3); % plot neck profile (curved)
                        plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:)+offset,'-k','LineWidth',1,'HandleVisibility','off') % right
                        hold on
                        plot(-MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:)+offset,'-k','LineWidth',1,'HandleVisibility','off') % left
                    elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                        % extracting data points for low curved necks
                        subplot(1,3,2)
                        plot(MS(n,m).Result(k,p).ExcessArea,MS(n,m).Result(k,p).DeltaSA,'.k','MarkerSize',13,'HandleVisibility','off') % mark point
                        hold on
                        sigma_text_flat{num_sigma_flat}=['$\tilde{\sigma}$ = ' num2str(MS(n,m).Result(k,p).sigmaD,2)];
                        num_sigma_flat=num_sigma_flat+1;
                        subplot(1,3,1); % plot neck profile (flat)
                        plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:)+offset,'-k','LineWidth',1,'HandleVisibility','off') % right
                        hold on
                        plot(-MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:)+offset,'-k','LineWidth',1,'HandleVisibility','off') % left
                    end
                end
            end
        end
    end
    offset=offset+delta;
end
%
subplot(1,3,1)
axis equal
xlim([-20 20])
ylim([-20 5])
text(0,0,sigma_text_flat,'Interpreter','latex')
%
subplot(1,3,3)
axis equal
xlim([-20 20])
ylim([-20 5])
text(0,0,sigma_text_curved,'Interpreter','latex')

%% plotting neck shape at cursor position in plot 1
evalin('base','N_click=1');
datacursormode on
dcm_obj = datacursormode(f1);
dcm_obj.updateDataCursors;
% Set update function
set(dcm_obj,'UpdateFcn',{@PlotNeckAtCursor_energy,MS,f1})
%
%% plotting (area,neck radius)
%% Figure 2
f2=figure();
subplot(1,2,1)
for n=1:N_tauD
    plot(ExcessAreaCurved{n},NeckRadiusCurved{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
    hold on
    flat_area_point(n)=plot(ExcessAreaFlat{n},NeckRadiusFlat{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
end
h=gca;

legend(flat_area_point,legend_text,'Position',[0.45 0.45 0.2 0.3]); % add new legend for each value of tauD
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5,'HandleVisibility','off');
%xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$')
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$','Interpreter','latex')
ylabel('Neck radius')
xlim([-70 50])
ylim([0 7])
A=linspace(0,abs(h.XLim(1)),10000);
r=sqrt(A/pi);
plot(-A,r,'--k','HandleVisibility','off');
%
%% plotting (excess area, r1)
f3=figure();
for n=1:N_tauD
    plot(ExcessAreaCurved{n},r1Curved{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
    hold on
    flat_area_point(n)=plot(ExcessAreaFlat{n},r1Flat{n},'.-','LineWidth',1,'Markersize',9,'Color',plot_color(n,:));
end
h=gca;
legend(flat_area_point,legend_text); % add new legend for each value of tauD
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5,'HandleVisibility','off');
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$','Interpreter','latex')
ylabel('edge radius r1','Interpreter','latex')

%% plotting neck shape at cursor position in plot 2
datacursormode on
dcm_obj = datacursormode(f2);
dcm_obj.updateDataCursors;
% Set update function
set(dcm_obj,'UpdateFcn',{@PlotNeckAtCursor_radius,MS,f2})