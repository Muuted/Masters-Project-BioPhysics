% This script produces a single solution for a given tauD and sigmaD
% value using shooting.
% Output is the Result structure containing the solution.
% Result is plotted in a separate script
%
clc
clearvars;
% input parameters:
alpha_target = -0.75;
tauD = 1;
sigmaD = 0.1;
%
% Calling ShapeAlpha to generate solution
[Result]=ShapeAlpha(alpha_target,tauD,sigmaD);
%
%
%% Plot all neck profiles within Result in separate subplots:
figure()
[R,C]=size(Result);
for k=1:R % loop over crossings t=tau
    for p=1:C % loop over psi2zero values
        if ~isempty(Result(k,p).ShapeSolution)
        subplot(R,C,C*(k-1)+p);
        plot(Result(k,p).Y(1,:),Result(k,p).Y(5,:)) % plot (r,z) in solution
        hold on;
        %plot(-Result(k,p).Y(1,:),Result(k,p).Y(5,:)) % plot mirrored (r,z) in solution
        plot(Result(k,p).ShapeSolution.ye(1,k),Result(k,p).ShapeSolution.ye(5,k),'og'); % plot integration endpoint
        axis equal
        xlabel('rD');
        ylabel('zD');
        plottext{1}=['\alpha_{out}=' num2str(Result(k,p).alpha_out)];
        plottext{2}=['r0D=' num2str(Result(k,p).r0D)];
        plottext{3}=['k=' num2str(k) '  p=' num2str(p) ];
        plottext{4}=['DeltaSA=' num2str(Result(k,p).DeltaSA)];
        plottext{5}=['Area=' num2str(Result(k,p).Area)];
        text(mean(xlim),mean(ylim)+0.25*range(ylim),plottext)
        
        m = length( Result(k,p).Y(1) );
        
        r_11 = Result(k,p).Y(1,m);
        psi_11 = Result(k,p).Y(2,m);
        dpsidt_11 = Result(k,p).Y(3,m);
        
        r_1 = Result(k,p).ShapeSolution.ye(1,k);
        psi_1 = Result(k,p).ShapeSolution.ye(2,k);
        dpsidt_1 = Result(k,p).ShapeSolution.ye(3,k);
        
        a11 = (1 - dpsidt_11)*r_11/sin(psi_11)  - 1 
        a1 = (1 - dpsidt_1)*r_1/sin(psi_1)  - 1 
        else
            continue
        end
    end
end





%% plot both neck shape solutions of Result in same plot:
figure()
[R,C]=size(Result);
for k=1:R % loop over crossings t=tau
    for p=1:C % loop over psi2zero values
        if ~isempty(Result(k,p).ShapeSolution)
        plot(Result(k,p).Y(1,:),Result(k,p).Y(5,:)) % plot (r,z) in solution
        hold on;
        %plot(-Result(k,p).Y(1,:),Result(k,p).Y(5,:)) % plot mirrored (r,z) in solution
        plot(Result(k,p).ShapeSolution.ye(1,k),Result(k,p).ShapeSolution.ye(5,k),'.g','MarkerSize',10); % plot integration endpoint
        axis equal
        xlabel('rD');
        ylabel('zD');
        plottext{1}=['\alpha_{out}=' num2str(Result(k,p).alpha_out)];
        plottext{2}=['r0D=' num2str(Result(k,p).r0D)];
        plottext{3}=['k=' num2str(k) '  p=' num2str(p) ];
%        text(mean(xlim),mean(ylim)+0.25*range(ylim),plottext)
        else
            continue
        end
    end
end
%
%
%% plot Hamiltonian in Result
% figure()
% [R,C]=size(Result);
% for k=1:R % loop over crossings t=tau
%     for p=1:C % loop over psi2zero values
%         if ~isempty(Result(k,p).ShapeSolution)
%         subplot(R,C,C*(k-1)+p);
%         plot(Result(k,p).SD,Result(k,p).HD,'-b'); % plot integration endpoint        
%         xlabel('SD');
%         ylabel('Hamiltonian');
%         legend(['\alpha_{out}=' num2str(Result(k,p).alpha_out)])
%         else
%             continue
%         end
%     end
% end
%
%
%
%% Old code to plot neck profiles in separate subplots. Using Y directly.
%
% set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
% set(groot, 'defaultAxesTickLabelInterpreter','Latex'); 
% set(groot, 'defaultLegendInterpreter','Latex');
% set(groot,'DefaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual');
% set(0,'defaultAxesFontSize',20);
% 
% figure()
% % set(groot,'defaulttextinterpreter','latex');
% % set(groot, 'defaultAxesTickLabelInterpreter','latex');  
% % set(groot, 'defaultLegendInterpreter','latex');
% % set(groot, 'defaultLegendFontSize',30);
% % set(groot, 'defaulttextFontSize',30);
% %set(groot, 'defaultAxesTickLabelFontSize',30);
% % set(groot,'defaulttextinterpreter','none');
% % set(groot, 'defaultAxesTickLabelInterpreter','none');  
% % set(groot, 'defaultLegendInterpreter','none');
% % 
% 
% %
% subplot(2,2,1)
% plot(Y(1,:),Y(5,:),'Linewidth',1) % plot (r,z) in solution
% title('Neck profile, one side')
% %xlabel('$$\tilde{r}$$');
% %ylabel('$$\tilde{z}$$');
% xlabel('Radius r');
% ylabel('Height z');
% axis equal
% %***********plot tension (s,t) *************
% subplot(2,2,2)
% plot(SD,Y(4,:),'m','Linewidth',2) % plot tension as function of s
% hold on
% plot(SD,HD,'r','Linewidth',2) % plot hamiltonian
% ylim([-0.1*max(Y(4,:)) inf])
% %xlabel('$$\tilde{s}$$');
% %ylabel('$$\tilde{t}$$');
% xlabel('Arc length $\tilde{s}$','Interpreter','Latex');
% ylabel('Tensions $\tilde{t}$, $\tilde{\mathcal{H}}$','Interpreter','latex');
% %title('Tension and Hamiltonian')
% legend('Tension, $\tilde{t}$','Hamiltonian $\tilde{\mathcal{H}}$','Location','east','Interpreter','latex');
% set(gca,'LineWidth',1);
% %legend('boxoff')
% %************plot (r,z) symmetric around hole ***********
% subplot(2,2,1)
% p1=plot(Y(1,:),Y(5,:),'k','Linewidth',2); % profile
% hold on
% xlim([-10 10]);
% ylim([-10 10]);
% axis equal
% p2=plot(-Y(1,:),Y(5,:),'k','Linewidth',2); % plot profile other side
% p3=plot(r0D,0,'.r','Markersize',8);
% p4=plot(-r0D,0,'.r','Markersize',8);
% %xlabel('$$\tilde{r}$$');
% %ylabel('$$\tilde{z}$$');
% xlabel('Radius $\tilde{r}$','Interpreter','Latex');
% ylabel('Height $\tilde{z}$','Interpreter','Latex');
% %title('Membrane profile')
% set(gca,'LineWidth',1);
% legend([p1 p3],'Membrane profile',['Hole radius=' num2str(round(r0D,1))],'Location','northeast','Interpreter','Latex')
% strT=['$\tau$ = ' num2str(tauD) newline ...
%       '$\sigma$ = ' num2str(sigmaD) newline ...
%       '$\Delta G$ = ' num2str(round(DeltaG,2)) newline...
%     '$\psi_1$ = ' num2str(round(Y(2,end)*180/pi,1)) '$^{\circ}$'];
% text(-10, 10,strT,'FontSize',11,'Interpreter','Latex');
% %legend('boxoff')
% %************* plot curvatures
% C1=Y(3,:); 
% C2=sin(Y(2,:))./Y(1,:);
% SumC=C1+C2;
% subplot(2,2,3)
% plot(SD,C1,'r','Linewidth',2); % profile
% hold on 
% plot(SD,C2,'b','Linewidth',2);
% plot(SD,SumC,'g','Linewidth',2)
% %xlabel('$$\tilde{s}$$');
% %ylabel('Curvature in units of $1/c_0$');
% xlabel('Arc length $\tilde{s}$','Interpreter','Latex');
% ylabel('Curvatures','Interpreter','Latex');
% %title('Curvatures')
% set(gca,'LineWidth',1);
% legend('$1/R_1$','$1/R_2$','$2H=1/R_1+1/R_2$','northeast')
% %******************* plot 3D
% %subplot(2,2,4)
% figure()
% S_plot=linspace(0.5*NeckSolutionD(m2).xe,NeckSolutionD(m2).xe,25); % arc vector for 3D plotting
% YP=deval(NeckSolutionD(m2),S_plot); % evaluate solution in points suitable for 3D plotting
% DrawNeck3D_v2(YP(1,:),YP(5,:),r0D_zero(m2)) % call 3D plot function
% ax = gca;               % get the current axis
% ax.Clipping = 'off'; % avoids clipping in 3D
% %zoom(1.8);