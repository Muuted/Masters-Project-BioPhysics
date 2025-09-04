%% script to plot the neck energy versus area
% clickable version that plots neck profile in separate subplot
%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(groot,'defaultAxesFontName','Helvetica')
set(groot,'DefaultLegendFontName','Helvetica')
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%
%% plot area,DeltaSA for fixed tauD. Only necks with all psi angles.
%
[N_tauD,N_sigmaD]=size(MS);
psi1_cut=-0.8*pi; % cutoff between large and small type necks
Area0=pi*(MS(1,1).Result(1,1).r2D)^2; % Area of flattened region with same radius and no hole or excess area
%n=1; % index to tauD
f1=figure();
subplot(1,2,1);
plot_color=jet(N_tauD);
for n=1:N_tauD
    for m=1:N_sigmaD % loop sigmaD
        [R,C]=size(MS(n,m).Result);
        if R~=0 && C~=0 % proceed if Results is not empty
            for k=1:R % loop over crossings t=tau
                for p=1:C % loop over psi2zero values
                    if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty
                        if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut
                            plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).DeltaSA,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor','none','LineWidth',2); % plot (AreaD,DeltaSA) endpoint curvature above psi1_cut
%                            plot(100*MS(n,m).Result(k,p).Area/Area0,MS(n,m).Result(k,p).DeltaSA,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor','none','LineWidth',2); % plot (AreaD,DeltaSA) endpoint curvature above psi1_cut
%                            plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).tauD,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor','none'); % plot (AreaD,DeltaSA) endpoint curvature above psi1_cut
                            hold on;
                        elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut
                            plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).DeltaSA,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor',plot_color(n,:)); % plot (AreaD,DeltaSA) endpoint curvature below psi1_cut
%                           plot(100*MS(n,m).Result(k,p).Area/Area0,MS(n,m).Result(k,p).DeltaSA,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor',plot_color(n,:)); % plot (AreaD,DeltaSA) endpoint curvature below psi1_cut
%                           plot(MS(n,m).Result(k,p).Area,MS(n,m).Result(k,p).tauD,'o-','MarkerSize',4,'MarkerEdgeColor',plot_color(n,:),'MarkerFaceColor',plot_color(n,:)); % plot (AreaD,DeltaSA) endpoint curvature below psi1_cut
                            hold on;
                        end
                    end
                end
            end
        end
    end
end
% Add vertical line at Area=zero hole radius
h=gca;
line([Area0 Area0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
%line([100 100],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
xlabel('Neck area (unit=c_0^{-2})')
%xlabel('Relative area (%)')
ylabel('Energy difference \DeltaS=S_{neck}-S_{flat}')
title('Neck energy versus Area')
%% plotting neck shape at cursor position in plot 1
datacursormode on
dcm_obj = datacursormode(f1);
dcm_obj.updateDataCursors;
% Set update function
set(dcm_obj,'UpdateFcn',{@PlotNeckAtCursor,MS,f1}) 
%% plot area,DeltaSA for fixed tauD. Only plot neck with largest psi angle
% procedure below is needed to sort empty array elements and select
% larges curvature
%
% [N_tauD,N_sigmaD]=size(MS);
% n=1; % index to tauD
% figure()
% for m=1:N_sigmaD % loop sigmaD
%     [R,C]=size(MS(n,m).Result);
%     if R~=0 && C~=0 % proceed if Results is not empty
%         psi1_largest=0;
%         for k=1:R % loop over crossings t=tau
%             for p=1:C % loop over psi2zero values
%                 if ~isempty(MS(n,m).Result(k,p).ShapeSolution) % proceed if ShapeSolution is not empty
%                     if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_largest % select most curved neck
%                         psi1_largest=MS(n,m).Result(k,p).ShapeSolution.ye(2,k);
%                         kc=k; % indices to must Result with most curved neck
%                         pc=p;
%                     end
%                 end
%             end
%         end
%         %        plot(MS(n,m).Result(kc,pc).sigmaD,MS(n,m).Result(kc,pc).Y(6,end),'.'); % plot (sigmaD,areaD)
%         %        plot(MS(n,m).Result(kc,pc).Y(6,end),MS(n,m).Result(kc,pc).SA,'.'); % plot (sigmaD,areaD)
%         plot(MS(n,m).Result(kc,pc).Area,MS(n,m).Result(kc,pc).DeltaSA,'.g'); % plot (AreaD,DeltaSA)
%         %plot(MS(n,m).Result(kc,pc).Y(1,:),MS(n,m).Result(kc,pc).Y(5,:));
%         hold on;
%         %        legendtext{m}=['$\sigma_D$ = ' num2str(MS(n,m).Result(kc,pc).sigmaD)];
%     else
%         continue
%     end
% end
% h=gca;
% %line([0 0],h.YLim,'Color','k','LineStyle',':','LineWidth',1)
% %axis equal
% %ylim(h.YLim);
% %xlim(h.XLim-1);
% %legend(legendtext);
% xlabel('Area')
% ylabel('DeltaSA')
