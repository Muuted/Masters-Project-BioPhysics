%% script to make a combined plot of content of Result using subplots
% Input to this script is the Result structure which must be generated in 
% for figure for paper
% 
% formatting figure
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');  
set(groot,'DefaultAxesTitleFontWeight','normal');
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%
%
figure()
% Input: select two solutions using two pairs of k and p values: (k1,p1)
% and (k2,p2)
% neck solution 1;
k1=2;
p1=1;
% neck solution 2:
% k2=4;
% p2=1;
k2=2;
p2=2;
%
%% plot shooting
subplot(2,2,1)
for v=1:size(Result(k1,p1).alpha_deviation,1)
    semilogx(Result(k1,p1).psi2_vector,Result(k1,p1).alpha_deviation(v,:),'.k') %plotting deviation
    hold on;
end
line(xlim,[0 0],'Color','k');
plot(Result(k1,p1).psi2_zeros,zeros(size(Result(k1,p1).psi2_zeros)),'o','MarkerFaceColor','r','MarkerEdgeColor','k');
ylim([-10 10]);
xlim([min(Result(k1,p1).psi2_vector) max(Result(k1,p1).psi2_vector)]);
%set(gca,'XTick',min(Result(k1,p1).psi2_vector):max(Result(k1,p1).psi2_vector)/10:max(Result(k1,p1).psi2_vector)); 
xlabel('Start angle of integration $\psi_2$ (rad)')
ylabel('Relative deviation ($\alpha$-$\alpha_{T}$/$\alpha_{T}$)')
title('Shooting of $\psi_2$ to obtain $\alpha_T=k_G/k=-0.75$')
hsp1 = get(gca, 'Position'); 
%% plotting s,t
subplot(2,2,2)
% solution no 1
plot(Result(k1,p1).ShapeSolution.x,Result(k1,p1).ShapeSolution.y(4,:),'--k'); % plot t parameter full range
hold on
plot(Result(k1,p1).SD,Result(k1,p1).Y(4,:),'r','LineWidth',2); % plot t parameter until alpha detected (physical solution)
%
plot(Result(k1,p1).ShapeSolution.xe(1,1:end),Result(k1,p1).ShapeSolution.ye(4,1:end),'ok'); % Mark all crossings t=tau
plot(Result(k1,p1).ShapeSolution.xe(1,k1),Result(k1,p1).ShapeSolution.ye(4,k1),'ok','MarkerFaceColor','k'); % Mark crossing with alpha_deviation=0
% solution no 2
plot(Result(k2,p2).ShapeSolution.x,Result(k2,p2).ShapeSolution.y(4,:),'--k'); % plot t parameter full range
hold on
plot(Result(k2,p2).SD,Result(k2,p2).Y(4,:),'b','LineWidth',2); % plot t parameter until alpha detected (physical solution)
plot(Result(k2,p2).ShapeSolution.xe(1,1:end),Result(k2,p2).ShapeSolution.ye(4,1:end),'ok') % mark all crossings
plot(Result(k2,p2).ShapeSolution.xe(1,k2),Result(k2,p2).ShapeSolution.ye(4,k2),'ok','MarkerFaceColor','k'); % Mark crossing with alpha_deviation=0
%
line(xlim,[Result(k1,p1).tauD Result(k1,p1).tauD])
xlabel('Arc length $\tilde{s}$ ');
ylabel('Lagrange function $\tilde{t}(\tilde{s})$');
title('Lagrange function $\tilde{t}$')
%% plotting neck shape rD,zD
subplot(2,2,3)
% solution no 1
plot(Result(k1,p1).Y(1,:),Result(k1,p1).Y(5,:),'r','LineWidth',2); % plot neck shape
hold on
plot(Result(k1,p1).ShapeSolution.ye(1,1:k1),Result(k1,p1).ShapeSolution.ye(5,1:k1),'ok') % all points with t=tau (open)
plot(Result(k1,p1).ShapeSolution.ye(1,k1),Result(k1,p1).ShapeSolution.ye(5,k1),'ok','MarkerFaceColor','k') % points with alpha=alpha_target (filled)
% solution no2
plot(Result(k2,p2).Y(1,:),Result(k2,p2).Y(5,:),'b','LineWidth',2); % plot neck shape
plot(Result(k2,p2).ShapeSolution.ye(1,1:k2),Result(k2,p2).ShapeSolution.ye(5,1:k2),'ok') % all points with t=tau (open)
plot(Result(k2,p2).ShapeSolution.ye(1,k2),Result(k2,p2).ShapeSolution.ye(5,k2),'ok','MarkerFaceColor','k') % points with alpha=alpha_target (filled)
%
xlabel('Radius $\tilde{r}$');
ylabel('Height $\tilde{z}$');
h=gca;
line([0 0],h.YLim,'Color','k','LineStyle','--','LineWidth',1)
axis equal
ylim([-0.2 2]);
xlim([-0.2 20]);
%xlim([-1 max(Result(k1,p1).Y(1,:))]);
plottext{1}=['$\alpha_{out}$=' num2str(Result(k1,p1).alpha_out)];
title('Neck profile')
set(gca, 'TitleFontWeight','normal');
%text(mean(xlim),0.7*mean(ylim),plottext)
hsp2 = get(gca, 'Position');                  
set(gca, 'Position', [hsp2(1:2)  hsp1(3) hsp2(4) ]) 
%% plot (s,mean curvature)
subplot(2,2,4)
hold on
% solution no 1
cp(1)=plot(Result(k1,p1).SD,sin(Result(k1,p1).Y(2,:))./Result(k1,p1).Y(1,:),'--r','LineWidth',1); % plot C1=1/R1
cp(2)=plot(Result(k1,p1).SD,Result(k1,p1).Y(3,:),':r','LineWidth',1); % plot C2=1/R2
cp(3)=plot(Result(k1,p1).SD,...
    (Result(k1,p1).Y(3,:) + sin(Result(k1,p1).Y(2,:))./Result(k1,p1).Y(1,:)),'r','LineWidth',2); % plot total curvature
plot(Result(k1,p1).ShapeSolution.xe(1,1:k1),...
    (Result(k1,p1).ShapeSolution.ye(3,1:k1) + sin(Result(k1,p1).ShapeSolution.ye(2,1:k1))./Result(k1,p1).ShapeSolution.ye(1,1:k1)),'ok'); % plot t=tau points
plot(Result(k1,p1).ShapeSolution.xe(1,k1),...
    (Result(k1,p1).ShapeSolution.ye(3,k1) + sin(Result(k1,p1).ShapeSolution.ye(2,k1))./Result(k1,p1).ShapeSolution.ye(1,k1)),'ok','MarkerFaceColor','k'); % plot endpoint
% solution no 2
cp(4)=plot(Result(k2,p2).SD,sin(Result(k2,p2).Y(2,:))./Result(k2,p2).Y(1,:),'--b','LineWidth',1); % plot C1=1/R1
cp(5)=plot(Result(k2,p2).SD,Result(k2,p2).Y(3,:),':b','LineWidth',1); % plot C2=1/R2
cp(6)=plot(Result(k2,p2).SD,...
    (Result(k2,p2).Y(3,:) + sin(Result(k2,p2).Y(2,:))./Result(k2,p2).Y(1,:)),'b','LineWidth',2); % total curvature
hold on
plot(Result(k2,p2).ShapeSolution.xe(1,1:k2),...
    (Result(k2,p2).ShapeSolution.ye(3,1:k2) + sin(Result(k2,p2).ShapeSolution.ye(2,1:k2))./Result(k2,p2).ShapeSolution.ye(1,1:k2)),'ok'); % t=tau points
plot(Result(k2,p2).ShapeSolution.xe(1,k2),...
    (Result(k2,p2).ShapeSolution.ye(3,k2) + sin(Result(k2,p2).ShapeSolution.ye(2,k2))./Result(k2,p2).ShapeSolution.ye(1,k2)),'ok','MarkerFaceColor','k'); % endpoint
% 
xlabel('Arc length $\tilde{s}$ ');
ylabel('Total curvature $\tilde{2H}=(\frac{1}{\tilde{R_1}}+\frac{1}{\tilde{R_2}}) $');
title('Mean curvature $\tilde{H}$');
%
legend_text_curvature{1}=['$C_1$'];
legend_text_curvature{2}=['$C_2$'];
legend_text_curvature{3}=['$C_1+C_2$'];
legend_text_curvature{4}=['$C_1$'];
legend_text_curvature{5}=['$C_2$'];
legend_text_curvature{6}=['$C_1+C2$'];
legend(cp,legend_text_curvature);
%
output_txt{1}=['sigmaD=' num2str(Result(k1,p1).sigmaD)];
output_txt{2}=['tauD=' num2str(Result(k1,p1).tauD)];
text(mean(xlim),mean(ylim),output_txt)
%% Plot profiles in 3D in separate figure windows
NP=25;
figure();
SurfY_3D(Result,k1,p1,NP)
% 
figure();
SurfY_3D(Result,k2,p2,NP)
%% Extra plots for supplementary
figure();
%% plot (s,curvature)
subplot(3,1,1)
hold on
box on
% solution no 1
cp(1)=plot(Result(k1,p1).SD,sin(Result(k1,p1).Y(2,:))./Result(k1,p1).Y(1,:),'--r','LineWidth',1); % plot C1=1/R1
cp(2)=plot(Result(k1,p1).SD,Result(k1,p1).Y(3,:),':r','LineWidth',1); % plot C2=1/R2
cp(3)=plot(Result(k1,p1).SD,...
    (Result(k1,p1).Y(3,:) + sin(Result(k1,p1).Y(2,:))./Result(k1,p1).Y(1,:)),'r','LineWidth',2); % plot total curvature
plot(Result(k1,p1).ShapeSolution.xe(1,1:k1),...
    (Result(k1,p1).ShapeSolution.ye(3,1:k1) + sin(Result(k1,p1).ShapeSolution.ye(2,1:k1))./Result(k1,p1).ShapeSolution.ye(1,1:k1)),'ok'); % plot t=tau points
plot(Result(k1,p1).ShapeSolution.xe(1,k1),...
    (Result(k1,p1).ShapeSolution.ye(3,k1) + sin(Result(k1,p1).ShapeSolution.ye(2,k1))./Result(k1,p1).ShapeSolution.ye(1,k1)),'ok','MarkerFaceColor','k'); % plot endpoint
% solution no 2
cp(4)=plot(Result(k2,p2).SD,sin(Result(k2,p2).Y(2,:))./Result(k2,p2).Y(1,:),'--b','LineWidth',1); % plot C1=1/R1
cp(5)=plot(Result(k2,p2).SD,Result(k2,p2).Y(3,:),':b','LineWidth',1); % plot C2=1/R2
cp(6)=plot(Result(k2,p2).SD,...
    (Result(k2,p2).Y(3,:) + sin(Result(k2,p2).Y(2,:))./Result(k2,p2).Y(1,:)),'b','LineWidth',2); % total curvature
hold on
plot(Result(k2,p2).ShapeSolution.xe(1,1:k2),...
    (Result(k2,p2).ShapeSolution.ye(3,1:k2) + sin(Result(k2,p2).ShapeSolution.ye(2,1:k2))./Result(k2,p2).ShapeSolution.ye(1,1:k2)),'ok'); % t=tau points
plot(Result(k2,p2).ShapeSolution.xe(1,k2),...
    (Result(k2,p2).ShapeSolution.ye(3,k2) + sin(Result(k2,p2).ShapeSolution.ye(2,k2))./Result(k2,p2).ShapeSolution.ye(1,k2)),'ok','MarkerFaceColor','k'); % endpoint
% 
xlabel('Arc length $\tilde{s}$ ');
ylabel('Curvature function');
%title('Curvature $\tilde{2H}$');
%
legend_text_curvature{1}=['$C_1$'];
legend_text_curvature{2}=['$C_2$'];
legend_text_curvature{3}=['$C_1+C_2$'];
legend_text_curvature{4}=['$C_1$'];
legend_text_curvature{5}=['$C_2$'];
legend_text_curvature{6}=['$C_1+C2$'];
legend(cp,legend_text_curvature);
%
output_txt{1}=['sigmaD=' num2str(Result(k1,p1).sigmaD)];
output_txt{2}=['tauD=' num2str(Result(k1,p1).tauD)];
text(mean(xlim),mean(ylim),output_txt)
%% plot s,psi
subplot(3,1,2)
box on
% solution no 1
plot(Result(k1,p1).SD,Result(k1,p1).Y(2,:),'r','LineWidth',2);
hold on
plot(Result(k1,p1).ShapeSolution.xe(1,1:k1),Result(k1,p1).ShapeSolution.ye(2,1:k1),'ok')
plot(Result(k1,p1).ShapeSolution.xe(1,k1),Result(k1,p1).ShapeSolution.ye(2,k1),'ok','MarkerFaceColor','k')
% solution no 2
plot(Result(k2,p2).SD,Result(k2,p2).Y(2,:),'b','LineWidth',2);
plot(Result(k2,p2).ShapeSolution.xe(1,1:k2),Result(k2,p2).ShapeSolution.ye(2,1:k2),'ok')
plot(Result(k2,p2).ShapeSolution.xe(1,k2),Result(k2,p2).ShapeSolution.ye(2,k2),'ok','MarkerFaceColor','k')
% 
xlabel('Arc length $\tilde{s}$');
ylabel('Neck angle $\psi$ (rad)');
%title('Membrane angle $\psi$');
%
output_txt{1}=['sigmaD=' num2str(Result(k1,p1).sigmaD)];
output_txt{2}=['tauD=' num2str(Result(k1,p1).tauD)];
text(mean(xlim),mean(ylim),output_txt)
%
%% plot Hamiltonian in Result
subplot(3,1,3)
box on
[R,C]=size(Result);
% plot solution 1
plot(Result(k1,p1).SD,Result(k1,p1).HD,'-r'); % plot integration endpoint  
hold on
% plot solution 2
plot(Result(k2,p2).SD,Result(k2,p2).HD,'-b'); % plot integration endpoint  
%
xlabel('Arc length $\tilde{s}$');
ylabel('Hamiltonian $\tilde{\mathcal{H}}$');
ylim([10*min(Result(k,p).HD) 10*max(Result(k,p).HD)]);
