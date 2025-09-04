% Generate first a solution to perturb:
clear all;
% input parameters
alpha_target=-0.75;
tauD=1;
sigmaD=-0.2;
%
[Result]=ShapeAlpha(alpha_target,tauD,sigmaD);
% plot the neck shape of the unperturbed solutions
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
            %xlim([0 max(Result(k,p).ShapeSolution.ye(1,k))]);
            xlabel('rD');
            ylabel('zD');
            plottext{1}=['\alpha_{out}=' num2str(Result(k,p).alpha_out)];
            plottext{2}=['r0D=' num2str(Result(k,p).r0D)];
            plottext{3}=['k=' num2str(k) '  p=' num2str(p) ];
            text(mean(xlim),0.7*mean(ylim),plottext)
        else
            continue
        end
    end
end
%%
%
NPoints=5000;
dSigma=0.1; % perturbation of tension SigmaD (to be added)
k=4;
p=1;
%
[SDp,Yp,SD,Y,alpha_p] = PerturbSigma(Result,k,p,NPoints,dSigma); % determining perturbed solution
%
[SAp] = Shape2Action(Yp,SDp,Result(k,p).tauD,Result(k,p).sigmaD,alpha_target); % energy of perturbed state (using original energy functional)
%[SAp] = Shape2Action(Yp,SDp,Result(k,p).tauD,Result(k,p).sigmaD+dSigma,alpha_target); % energy of perturbed state
[SA] = Shape2Action(Y,SD,Result(k,p).tauD,Result(k,p).sigmaD,alpha_target);
DSAf_sigma=(SAp-SA)/SA;
%
figure()
plot(Y(1,:),Y(5,:),'k','LineWidth',2) % plot (r,z) in original solution
hold on
plot(Yp(1,:),Yp(5,:),'b','LineWidth',2); % plot (r,z) in perturbed solution
title('Original and perturbed solution')
xlabel('rD');
ylabel('zD');
axis equal;
%
% plot areas
figure()
plot(SD,abs(Y(6,:)),'.-k')
hold on
plot(SDp,abs(Yp(6,:)),'.-b')
% plot psip
figure()
plot(SD,Y(2,:),'k')
hold on
plot(SDp,Yp(2,:),'b')
%plot(psipdot)
figure()
plot(SD,Y(3,:),'k')
hold on
plot(SDp,Yp(3,:),'.-b')
%