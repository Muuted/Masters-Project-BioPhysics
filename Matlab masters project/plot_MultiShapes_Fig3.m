%% script to plot the tau,sigma variation of shapes. Plots larghe and small
% curvature solutions in separate plots
% Used for figure 3 in paper
% - need as input the MS structure made with the MultiShapes_script.
%
% setting default options
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex')
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(groot,'DefaultLegendFontSize',16,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',25);
%
psi1_cut=-pi; % cutoff between large and small type necks
figure()
rp_max=5; % max radius for plotting
zp_max=2.5;
%% plot sigmaD variation for fixed tauD
MS=MS_sigma;
n=1; % index to tauD
[N_tauD,N_sigmaD]=size(MS);
plot_color=jet(N_sigmaD); % color of points
for m=1:N_sigmaD % loop sigmaD
    [R,C]=size(MS(n,m).Result);
    if R~=0 && C~=0
        psi1_largest=0;
        for k=1:R % loop over crossings t=tau
            for p=1:C % loop over psi2zero values
                if ~isempty(MS(n,m).Result(k,p).ShapeSolution)
                    if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut % select most curved neck
                        subplot(2,2,1)
                        plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:),'LineWidth',2,'Color',plot_color(m,:));
                        hold on;
                        lg1_text{m}=['$\tilde{\sigma}$ = ' num2str(MS(n,m).Result(k,p).sigmaD,2) ', Excess area = ' num2str(MS(n,m).Result(k,p).ExcessArea,2)];
                        subplot(2,2,[2,4])
                        hold on;
                        ExcessAreaCurved(m)=MS(n,m).Result(k,p).ExcessArea;
                        DeltaSACurved(m)=MS(n,m).Result(k,p).DeltaSA;
                        plot(MS(n,m).Result(k,p).ExcessArea,MS(n,m).Result(k,p).DeltaSA,'^-','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',plot_color(m,:));
                    elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut %
                        subplot(2,2,3)
                        plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:),'LineWidth',2,'Color',plot_color(m,:));
                        hold on;                                
                        lg3_text{m}=['$\tilde{\sigma}$ = ' num2str(MS(n,m).Result(k,p).sigmaD,2) ', Excess area = ' num2str(MS(n,m).Result(k,p).ExcessArea,2)];
                        subplot(2,2,[2,4])
                        hold on;
                        plot(MS(n,m).Result(k,p).ExcessArea,MS(n,m).Result(k,p).DeltaSA,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',plot_color(m,:));
                    end
                end
            end
        end
    else
        continue
    end
end
%
subplot(2,2,[2,4])
hold on
plot(ExcessAreaCurved,DeltaSACurved,'-k');
%
subplot(2,2,1)
h=gca;
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1)
axis equal
ylim([-0.2 zp_max]);
xlim([-0.2 rp_max]);
legend(lg1_text);
xlabel('Radius $\tilde{r}$')
ylabel('Height $\tilde{z}$')
title('Curved necks ($\psi_1 < -\pi$) with $\tilde{\tau}$=1')
box on;
ax=gca;
ax.LineWidth=2;
%
subplot(2,2,3)
h=gca;
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1)
axis equal
ylim([-0.2 zp_max]);
xlim([-0.2 rp_max]);
xlabel('Radius $\tilde{r}$')
ylabel('Height $\tilde{z}$')
legend(lg3_text);
title('Flat necks ($\psi_1 > -\pi$) with $\tilde{\tau}$=1')
box on;
ax=gca;
ax.LineWidth=2;
%
subplot(2,2,[2,4])
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1.5);
xlabel('Excess area $\Delta \tilde{A} = \tilde{A}_{neck}-\tilde{A}_{disc}$')
ylabel('Relative energy  $\Delta \tilde{S}=\tilde{S}_{neck}-\tilde{S}_{flat}$')
title('Relative neck energy versus excess area ($\tilde{\tau}$=1)')
box on;
ax=gca;
ax.LineWidth=2;
%
%% plot tauD variation for fixed sigmaD
% MS=MS_tau;
% m=1; % index to fixed sigmaD
% [N_tauD,N_sigmaD]=size(MS);
% %figure()
% for n=1:N_tauD % loop tauD
%     [R,C]=size(MS(n,m).Result);
%     if R~=0 && C~=0
%         psi1_largest=0;
%         for k=1:R % loop over crossings t=tau
%             for p=1:C % loop over psi2zero values
%                 if ~isempty(MS(n,m).Result(k,p).ShapeSolution)
%                     if MS(n,m).Result(k,p).ShapeSolution.ye(2,k)<psi1_cut % select most curved neck
%                         subplot(2,2,2)
%                         plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:),'LineWidth',2);
%                         hold on;
%                         lg2_text{n}=['$\tilde{\tau}$ = ' num2str(MS(n,m).Result(k,p).tauD,2) ', Excess area = ' num2str(MS(n,m).Result(k,p).ExcessArea,2)];
%                     elseif MS(n,m).Result(k,p).ShapeSolution.ye(2,k)>psi1_cut %
%                         subplot(2,2,4)
%                         plot(MS(n,m).Result(k,p).Y(1,:),MS(n,m).Result(k,p).Y(5,:),'LineWidth',2);
%                         hold on;                                
%                         lg4_text{n}=['$\tilde{\tau}$ = ' num2str(MS(n,m).Result(k,p).tauD,2) ', Excess area = ' num2str(MS(n,m).Result(k,p).ExcessArea,2)];
%                     end
%                 end
%             end
%         end
%     else
%         continue
%     end
% end
% %
% subplot(2,2,2)
% h=gca;
% line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1)
% axis equal
% ylim([-0.2 zp_max]);
% xlim([-0.2 rp_max]);
% legend(lg2_text);
% xlabel('Radius $\tilde{r}$')
% ylabel('Height $\tilde{z}$')
% title('Solution 1 ($\psi_1 < -\pi$)')
% %
% subplot(2,2,4)
% h=gca;
% line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1)
% axis equal
% ylim([-0.2 zp_max]);
% xlim([-0.2 rp_max]);
% xlabel('Radius $\tilde{r}$')
% ylabel('Height $\tilde{z}$')
% legend(lg4_text);
% title('Solution 2 ($\psi_1 > -\pi$)')

