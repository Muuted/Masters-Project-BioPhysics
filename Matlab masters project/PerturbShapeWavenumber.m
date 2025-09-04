% Function to perturb the neck shape using a sinusoidal perturbation with
% defined amplitude (b) wavenumber (wn) and phase (ph)
%
% b=amplitude of perturbations in units= max(z)
% wn=wavenumber of perturbations in units=1/range(SD)
% ph=phase angle of sin perturbation (unit=rad)
%
function [SDp,Yp,SD,Y,alpha_p] = PerturbShapeWavenumber(ShapeSolution,k,NPoints,b,wn,ph)
% First evaluate the pathlength endpoint fe that gives same area of
% perturbed and unperturbed solutions
zerofun = @(fe) ExtendSolution(ShapeSolution,k,fe,NPoints,b,wn,ph); % function handle. Evaluates area difference between perturbed and original solution as function of fe (extension parameter)
% Note: argument k is the index (enumeration) of the crossing point t=tau. Is needed when evaluating and
% picking the correct solution from ShapeSolution which contains all the crossings.
%
fe0=fzero(zerofun,[0.4 1.4]); % evaluates the zero point value of fe
%
% After fe0 is found we evaluate the final perturbed solution for output
[~,rp,zp,psip,psipdot,SD,Y,SD_tot,~,area_pert] = ExtendSolution(ShapeSolution,k,fe0,NPoints,b,wn,ph); % call ExtendRange with solution fe0
% placing perturbed functions in the Yp matrix (same format as original solution)
Yp(1,:)=rp;
Yp(2,:)=psip;
Yp(3,:)=psipdot;
Yp(4,:)=NaN; % we did not calculate t for the perturbed solution (not needed)
Yp(5,:)=zp;
Yp(6,:)=area_pert;
SDp=SD_tot;
%
% find alpha in perturbed state (will not have the target value!,- does it matter?)
r1D=Yp(1,end); % radius at free edge
psi1=Yp(2,end); % angle at free edge
psidot1D=Yp(3,end); % dpsi/ds at free edge
alpha_p=(1-psidot1D)*r1D/sin(psi1)-1; % alpha for perturbed solution calculated as for the ordinary solution
%
% Start defining nested utility functions.
    function [delta_area,rp,zp,psip,psipdot,SD,Y,SD_tot,area_unpert,area_pert] = ExtendSolution(ShapeSolution,k,fe,NPoints,b,wn,ph)
        % ExtendSolution is total function to generate extended solution with perturbation.
        % ExtendSolution is called by the minimization function zerofun to achieve equal area relative to unperturbed solution.
        %
        if fe>1 % if longer path
            SD=linspace(ShapeSolution.x(1),ShapeSolution.xe(k),NPoints); % unextended pathvector
            Y=deval(ShapeSolution,SD); % original unextended solution
            %
            SD_tot=linspace(ShapeSolution.x(1),fe*ShapeSolution.xe(k),NPoints); % extended pathvector
            Y_tot=deval(ShapeSolution,SD_tot); % extended solution (will produce error if SD_tot is beyond range in ShapeSolution)
            [rp,zp,psip,psipdot]=unpert2pert(Y_tot,SD_tot,b,wn,ph); % perturb extended solution
            area_unpert=rz2area(Y(1,:),Y(5,:)); % area unperturbed solution (vector)
            area_pert=rz2area(rp,zp); % area perturbed (vector)
            delta_area=area_pert(end)-area_unpert(end); % area difference (number). Must be minimized to zero by adjusting fe
        else % if fe<1 i.e. short path
            SD=linspace(ShapeSolution.x(1),ShapeSolution.xe(k),NPoints); % unextended pathvector
            Y=deval(ShapeSolution,SD); % original unextended solution
            %
            SD_tot=linspace(ShapeSolution.x(1),fe*ShapeSolution.xe(k),NPoints); % shortened pathvector
            Y_tot=deval(ShapeSolution,SD_tot); % evaluate solution over shortened range
            [rp,zp,psip,psipdot]=unpert2pert(Y_tot,SD_tot,b,wn,ph); % obtain perturbed solution
            area_unpert=rz2area(Y(1,:),Y(5,:)); % area unperturbed solution (vector)
            area_pert=rz2area(rp,zp); % area perturbed (vector)
            delta_area=area_pert(end)-area_unpert(end); % area difference (number). Must be minimized to zero by adjusting fe
        end
    end
%
% Second utility function: Definition of rz2area
    function [area]=rz2area(r,z)
        % Utility function to determine area of neck from r and z coordinates
        ds=sqrt(diff(z).^2 + diff(r).^2); % pathlength between points on neck profile
        dA=2*pi*ds.*0.5.*(r(2:end)+r(1:end-1)); % area increment betwen points. Using average radius between endpoints
        dA=[0 dA];
        area=cumsum(dA); % area as cummulative sum
    end
% 
% Third utility function: Definition of unpert2pert
    function [rp,zp,psip,psipdot]=unpert2pert(Y,SD,b,wn,ph)
        % utility function to generate perturbed neck shape (rp,zp) from unperturbed solution (Y).
        % SD = pathlength vector for solution Y
        % b=amplitude of perturbations in units= max(z)
        % l=wavenumber of perturbations 
        % ph=phase angle of sin perturbation (unit=rad)
        a=b*max(Y(5,:))*abs((SD/range(SD)).^3).*sin(((2*pi*wn)/(range(SD)))*SD + ph); % perturbation as function of SD
        zp=Y(5,:)+a.*cos(Y(2,:)); % perturbed z-coordinate
        rp=Y(1,:)-a.*sin(Y(2,:)); % perturbed r-coordinate
        % obtain psip = perturbed psi-coordinate
        dsp=sqrt(diff(zp).^2 + diff(rp).^2); % pathlength step between points on neck profile (positive)
        psip=atan(diff(zp)./diff(rp)); % must use tan to be able to unwrap angles !
        psip=0.5*unwrap(2*psip);
        %extrapolaring first point (missing because of diff):
        sp_startpoints=-cumsum(dsp(1:10)); % first 10 points in s
        psip2=interp1(sp_startpoints,psip(1:10),0,'linear','extrap'); % extrapolating to psi2
        psip=[psip2 psip]; % appending psi2 to psi
        dsp=[dsp(1) dsp]; % appending ds(1) to s
        % finding psipdot=dpsip/dsp or perturbed psidot
        psipdot=-diff(psip)./dsp(2:end); % calculating psidot by differentiation
        psidot2=interp1(sp_startpoints,psipdot(1:10),0,'linear','extrap'); % extrapolating to psi2
        psipdot=[psidot2 psipdot]; %appending psidot2 to psidot
    end
end