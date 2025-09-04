function output_txt = PlotNeckAtCursor(~,event_obj,MS,f1)
% ~            Currently not used (empty)
% event_obj    Object containing event data structure
% output_txt   Data cursor text
pos = get(event_obj, 'Position'); % get position of mouse click
disp(pos);
% start identifying the m,n,k,p value corresponding to the currently clicked point
[N_tauD,N_sigmaD]=size(MS);
%disp(N_tauD);disp(N_sigmaD)
for n=1:N_tauD
    for m=1:N_sigmaD % loop sigmaD
        [R,C]=size(MS(n,m).Result);
        for k=1:R % loop over crossings t=tau
            for p=1:C % loop over psi2zero values
                 if isequal(MS(n,m).Result(k,p).ExcessArea,pos(1)) && isequal(MS(n,m).Result(k,p).DeltaSA,pos(2))
                     mp=m; % idenitifying indices m,n,k,p to clicked point in plot
                     np=n;
                     kp=k;
                     pp=p;
                 end
            end
        end
    end
end
% stop identifying the m,n,k,p value corresponding to the currently clicked point
% making output text for plot
output_txt{1}=['ExcessArea=' num2str(MS(np,mp).Result(kp,pp).ExcessArea)];
output_txt{2}=['DeltaSA=' num2str(MS(np,mp).Result(kp,pp).DeltaSA)];
output_txt{3}=['sigmaD=' num2str(MS(np,mp).Result(kp,pp).sigmaD)];
output_txt{4}=['tauD=' num2str(MS(np,mp).Result(kp,pp).tauD)];
%
%% plot neck profile associated with a clicked point in the (Area,DeltaSA) plot,- using mp,np,kp,pp, identified above
figure(f1); % set to main figure
%
subplot(1,2,2); 
cla;
%rp_max=5;
%cla; % clear subplot axes
plot(MS(np,mp).Result(kp,pp).Y(1,:),MS(np,mp).Result(kp,pp).Y(5,:),'k','Linewidth',1); % plot profile right side 
hold on
plot(-MS(np,mp).Result(kp,pp).Y(1,:),MS(np,mp).Result(kp,pp).Y(5,:),'k','Linewidth',1); % plot profile left side
line([0 0],[0 3],'Color','k','LineStyle','--','LineWidth',1.5);
%xlim([-rp_max rp_max]);
axis equal
xlabel('Radius r');
ylabel('Height z');
title('Membrane profile')
YL=ylim;
text(mean(xlim),0.8*YL(1)+0.2*YL(2),output_txt);
