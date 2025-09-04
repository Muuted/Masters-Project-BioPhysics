% This script produces a single solution for a given tauD and sigmaD
% value using shooting.
% Output is the Result structure containing the solution.
% Result is plotted in a separate script
%
clear all;
clc
clearvars;
% input parameters:
alpha_target=-0.75;
tauD=1;
sigmaD=-0.2;
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
        else
            continue
        end
    end
end
%
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
