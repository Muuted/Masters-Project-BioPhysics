%% Script to plot the energy difference DeltaSA for the neck versus the neck area.
% Input is MS - structure with necks for varying tauD and sigmaD.
% MS is produced with MultiShapes_script 
%
% clickable version that plots neck profile in separate subplot when
% clicking on a point in the (Area,DeltaSA) plot.
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
%% plot (Area,DeltaSA) for fixed tauD. Necks with all psi angles are plotted.
%
[N_tauD,N_sigmaD]=size(MS);
psi1_cut=-pi; % cutoff between large and small type necks
f1=figure();
subplot(1,2,1);
plot_color=jet(N_tauD);
%plot_color='k';
for n=1:N_tauD
    for m=1:N_sigmaD % loop sigmaD
        [R,C]=size(MS(n,m).Result);
        if R~=0 && C~=0 % proceed if Result is not empty
            for k=1:R % loop over crossings t=tau
                for p=1:C % loop over psi2zero values
                    if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty
                        Area0=pi*(MS(n,m).Result(k,p).r2D)^2; % Area of flattened region with same radius and no hole or excess area
                        if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut
                            hpoint(n)=plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).DeltaSA,'o','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor','none','LineWidth',2); % plot (AreaD,DeltaSA) endpoint curvature above psi1_cut
                            hold on;
                        elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                            hpoint(n)=plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).DeltaSA,'o','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor',plot_color(n,:)); % plot (AreaD,DeltaSA) endpoint curvature below psi1_cut
                            hold on;
                            legend_text{n}=['$\tilde{\tau}$ = ' num2str(MS(n,m).Result(k,p).tauD,2)];
                        end
                    end
                end
            end
        end
    end
end
%
h=gca;
%line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
line([Area0 Area0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
legend(hpoint,legend_text);
%line([100 100],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
xticks([Area0-5 Area0 Area0+5 Area0+10 Area0+15 Area0+20 Area0+25 Area0+30 Area0+35 Area0+40])
xticklabels({'-5','0','5','10','15','20','25','30','35','40'})
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$')
ylabel('Relative energy  $\Delta \tilde{S}=\tilde{S}_{neck}-\tilde{S}_{flat}$')
title('Neck energy versus area (constant $\tilde{\tau}$=1)')
%% plotting neck shape at cursor position in plot 1
datacursormode on
dcm_obj = datacursormode(f1);
dcm_obj.updateDataCursors;
% Set update function
set(dcm_obj,'UpdateFcn',{@PlotNeckAtCursor,MS,f1}) 