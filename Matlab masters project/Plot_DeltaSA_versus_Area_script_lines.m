%% Script to plot the energy difference DeltaSA for the neck versus the neck area.
%
% Setting the default properties of the axis
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}) % set to real black
set(groot,'defaultAxesFontName','Helvetica')
set(groot,'DefaultLegendFontName','Helvetica')
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%
%% plot (Area,DeltaSA) for fixed tauD. 
%
[N_tauD,N_sigmaD]=size(MS);
psi1_cut=-pi; % cutoff between large and small type necks
f1=figure();
subplot(1,2,1);
plot_color=jet(N_tauD);
%plot_color='b';
AH=cell(1,N_tauD);
DSH=cell(1,N_tauD);
AL=cell(1,N_tauD);
DSL=cell(1,N_tauD);
Neck_radius=cell(1,N_tauD);
for n=1:N_tauD
    ah=1; % index to high curvature points
    al=1; % index to low curvature points
    for m=1:N_sigmaD % loop sigmaD
        [R,C]=size(MS(n,m).Result);
        if R~=0 && C~=0 % proceed if Result is not empty
            for k=1:R % loop over crossings t=tau
                for p=1:C % loop over psi2zero values
                    if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty
                        Area0=pi*(MS(n,m).Result(k,p).r2D)^2; % Area of flattened region with same radius and no hole or excess area
                        if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut
                            % extracting data points for high curved necks
                            AH{n}(ah)=MS(n,m).Result(k,p).Area; % neck area
                            DSH{n}(ah)=MS(n,m).Result(k,p).DeltaSA; % relative energy
                            Neck_radius{n}(ah)=min(MS(n,m).Result(k,p).Y(1,:)); % neck radius as minimum radius of shape
                            ah=ah+1;
                        elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                            % extracting data points for low curved necks
                            AL{n}(al)=MS(n,m).Result(k,p).Area;
                            DSL{n}(al)=MS(n,m).Result(k,p).DeltaSA;
                            al=al+1;
                            legend_text{n}=['$\tilde{\tau}$ = ' num2str(MS(n,m).Result(k,p).tauD,2)]; % construct legend text
                        end
                    end
                end
            end
        end
    end
end
% plotting extracted data points
for n=1:N_tauD
    plot(AH{n},DSH{n},'.-','LineWidth',2,'Color',plot_color(n,:));
    hold on
    hpoint(n)=plot(AL{n},DSL{n},'.-','LineWidth',2,'Color',plot_color(n,:));
end
% Add vertical line at Area=zero hole radius
h=gca;
line([Area0 Area0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
line([Area0+5 Area0+5],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
line([Area0+20 Area0+20],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
legend(hpoint,legend_text); % add new legend for each value of tauD
xticks([Area0-5 Area0 Area0+5 Area0+10 Area0+15 Area0+20 Area0+25 Area0+30 Area0+35 Area0+40])
xticklabels({'-5','0','5','10','15','20','25','30','35','40'})
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$')
ylabel('Relative energy  $\Delta \tilde{S}=\tilde{S}_{neck}-\tilde{S}_{flat}$')
title('Neck energy versus area')
%% plotting (area,neck radius)
figure()
for n=1:N_tauD
    plot(AH{n},Neck_radius{n},'-','LineWidth',2,'Color',plot_color(n,:));
    hold on
    % extrapolate neck radii to excess area=0
    % p=polyfit(AH{n}(1:10),Neck_radius{n}(1:10),1);
    % x=linspace(Area0,AH{n}(10),10);
    % y=polyval(p,x);
    % plot(x,y,'--k');
end
h=gca;
line([Area0 Area0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
line(h.XLim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5);
xticks([Area0-5 Area0 Area0+5 Area0+10 Area0+15 Area0+20 Area0+25 Area0+30 Area0+35 Area0+40])
xticklabels({'-5','0','5','10','15','20','25','30','35','40'})
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$')
ylabel('Neck radius')
 



%% plotting neck shape at cursor position in plot 1
datacursormode on
dcm_obj = datacursormode(f1);
dcm_obj.updateDataCursors;
% Set update function
set(dcm_obj,'UpdateFcn',{@PlotNeckAtCursor,MS,f1})