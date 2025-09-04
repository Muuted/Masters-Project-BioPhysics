% plot Result
%
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(groot, 'defaultLegendInterpreter','Latex');
set(groot,'DefaultAxesTitleFontWeight','normal');
set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
set(0,'defaultAxesFontSize',20);
%
figure()
k1=4;
p1=1;
k2=2;
p2=1;

%% plotting r,z
subplot(2,2,3)
plot(Result(k1,p1).Y(1,:),Result(k1,p1).Y(5,:),'g','LineWidth',2); % plot neck shape
axis equal
hold on
plot(Result(k1,p1).ShapeSolution.ye(1,1:k1),Result(k1,p1).ShapeSolution.ye(5,1:k1),'ok') % all points with t=tau (open)
plot(Result(k1,p1).ShapeSolution.ye(1,k1),Result(k1,p1).ShapeSolution.ye(5,k1),'ok','MarkerFaceColor','k') % point with alpha=alpha_target (filled)
% solution no2
plot(Result(k2,p2).Y(1,:),Result(k2,p2).Y(5,:),'g','LineWidth',2); % plot neck shape
plot(Result(k2,p2).ShapeSolution.ye(1,1:k2),Result(k2,p2).ShapeSolution.ye(5,1:k2),'ok') % all points with t=tau (open)
plot(Result(k2,p2).ShapeSolution.ye(1,k2),Result(k2,p2).ShapeSolution.ye(5,k2),'ok','MarkerFaceColor','k') % point with alpha=alpha_target (filled)
%
xlabel('r_D');
ylabel('z_D');
h=gca;
line([0 0],h.YLim,'Color','k','LineStyle',':','LineWidth',2)
axis equal
ylim(h.YLim);
xlim(h.XLim-1);
plottext{1}=['\alpha_{out}=' num2str(Result(k1,p1).alpha_out)];
title('Neck shape')
set(gca, 'TitleFontWeight','normal');
%text(mean(xlim),0.7*mean(ylim),plottext)
%% plotting s,t
subplot(2,2,2)
plot(Result(k1,p1).SD,Result(k1,p1).Y(4,:),'r','LineWidth',2); % plot t parameter
hold on
plot(Result(k1,p1).ShapeSolution.xe(1,1:k1),Result(k1,p1).ShapeSolution.ye(4,1:k1),'ok')
plot(Result(k1,p1).ShapeSolution.xe(1,k1),Result(k1,p1).ShapeSolution.ye(4,k1),'ok','MarkerFaceColor','k');
% solution no 2
plot(Result(k2,p2).SD,Result(k2,p2).Y(4,:),'r','LineWidth',2); % plot t parameter
plot(Result(k2,p2).ShapeSolution.xe(1,1:k2),Result(k2,p2).ShapeSolution.ye(4,1:k2),'ok')
plot(Result(k2,p2).ShapeSolution.xe(1,k2),Result(k2,p2).ShapeSolution.ye(4,k2),'ok','MarkerFaceColor','k');
%
line([minmax(Result(k1,p1).SD)],[Result(k1,p1).tauD Result(k1,p1).tauD])
xlabel('S_D');
ylabel('t_D');
title('Lagrange parameter t')
%% plot s,psi
subplot(2,2,4)
plot(Result(k1,p1).SD,Result(k1,p1).Y(2,:),'b','LineWidth',2);
hold on
plot(Result(k1,p1).ShapeSolution.xe(1,1:k1),Result(k1,p1).ShapeSolution.ye(2,1:k1),'ok')
plot(Result(k1,p1).ShapeSolution.xe(1,k1),Result(k1,p1).ShapeSolution.ye(2,k1),'ok','MarkerFaceColor','k')
% solution no 2
plot(Result(k2,p2).SD,Result(k2,p2).Y(2,:),'b','LineWidth',2);
plot(Result(k2,p2).ShapeSolution.xe(1,1:k2),Result(k2,p2).ShapeSolution.ye(2,1:k2),'ok')
plot(Result(k2,p2).ShapeSolution.xe(1,k2),Result(k2,p2).ShapeSolution.ye(2,k2),'ok','MarkerFaceColor','k')
% 
xlabel('S_D');
ylabel('\psi (rad)');
title('Membrane angle \psi')
%% plot shooting
subplot(2,2,1)
for v=1:size(Result(k1,p1).alpha_deviation,1)
    semilogx(Result(k1,p1).psi2_vector,Result(k1,p1).alpha_deviation(v,:),'.k') %plotting deviation
    hold on;
end
line(xlim,[0 0],'Color','k');
plot(Result(k1,p1).psi2_zeros,zeros(size(Result(k1,p1).psi2_zeros)),'o','MarkerFaceColor','r','MarkerEdgeColor','k');
xlabel('start angle \psi_2 (rad)')
ylabel('Relative deviation (\alpha-\alpha_{T}/\alpha_{T})')
title('Shooting to obtain \alpha=k_G/k=\alpha_T=-0.75')