% Function to perturb the neck shape using a sinusoidal perturbation with
% defined amplitude (b) wavelength (l) and phase (ph)
%
% b=amplitude of perturbations in units= max(z)
% l=wavelength of perturbations in units=range(SD)
% ph=phase angle of sin perturbation (unit=rad)
%
function [SDp,Yp,SD,Y,alpha_p] = PerturbSigma(Result,k,p,NPoints,dSigma)
% First evaluate the pathlength endpoint fe that gives same area of
% perturbed and unperturbed solutions
zerofun = @(fe) PerturbSigmaExtend(Result,k,p,fe,NPoints,dSigma); % function handle. Evaluates area difference between perturbed and original solution as function of fe (extension parameter)
% Note: argument k is the index (enumeration) of the crossing point t=tau. Is needed when evaluating and
% picking the correct solution from ShapeSolution which contains all the crossings.
%
fe0=fzero(zerofun,[0.4 1.4]); % evaluates the zero point value of fe
%
% After fe0 is found we evaluate the final perturbed solution for output
[~,SD,Y,SDp,Yp,alpha_p] = PerturbSigmaExtend(Result,k,p,fe0,NPoints,dSigma);
%
% Start defining nested utility functions.
    function [delta_area,SD,Y,SDp,Yp,alpha_p] = PerturbSigmaExtend(Result,k,p,fe,NPoints,dSigma)
        % PerturbSigmaExtend is total function to generate extended solution with perturbation of sigma.
        % PerturbSigmaExtend is called by the minimization function zerofun to achieve equal area relative to unperturbed solution.
        %
        [ShapeSolution_o,~,~]=Shape(Result(k,p).r2D,Result(k,p).psi2,Result(k,p).tauD,Result(k,p).sigmaD);
        % evaluate perturbed solution
        [ShapeSolution_p,alpha_p,~]=Shape(Result(k,p).r2D,Result(k,p).psi2,Result(k,p).tauD,Result(k,p).sigmaD + dSigma);
        %
        SD=linspace(ShapeSolution_o.x(1),ShapeSolution_o.xe(k),NPoints); % unextended pathvector
        Y=deval(ShapeSolution_o,SD); % original unextended solution
        %
        SDp=linspace(ShapeSolution_p.x(1),fe*ShapeSolution_p.xe(k),NPoints); % extended pathvector
        Yp=deval(ShapeSolution_p,SDp); % extended solution (will produce error if SD_tot is beyond range in ShapeSolution)
        %
        delta_area=Yp(6,end)-Y(6,end); % area difference (number). Must be minimized to zero by adjusting fe
    end
end