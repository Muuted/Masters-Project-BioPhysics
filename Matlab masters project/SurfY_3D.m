function SurfY_3D(Result,k,p,NP)
%%
rmaxp=7; % max radius for plotting
S_plot=linspace(rmaxp-Result(k,p).r2D,Result(k,p).ShapeSolution.xe(k),NP); % arc vector for 3D plotting
YP=deval(Result(k,p).ShapeSolution,S_plot); % evaluate solution in points suitable for 3D plotting
r_edge=YP(1,:);
z_edge=YP(5,:);
edge_angle=YP(2,:);
%% Creating membrane profile in 3D (X,Y,Z)
N_asm=16; % number of angles for azimuth rotation
phi=(0:pi/N_asm:2*pi)'; % angle for azimuth rotation
% Z is just a replication with the number of phi angles:
Z=z_edge(ones(1,length(phi)),:);
X=cos(phi)*r_edge(1,:);% X matrix is created by projecting r values using cos(phi)
Y=sin(phi)*r_edge(1,:);% Y matrix is created by projecting r values using sin(phi)
C=edge_angle(ones(1,length(phi)),:);
%% plotting of the curved membrane surface
% wih colored surface
% h=surf(X,Y,Z,C);
% colormap turbo;
% colorbar;
% with grey surface:
h=surf(X,Y,Z);
%set(h,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',1);
%set(h,'LineWidth',0.2);
set(h,'FaceColor',[1 1 1],'FaceAlpha',1);
set(h,'LineWidth',2);
hold on;
axis equal;
rotate3d on; % allow interactive rotation with mouse
axis off;
ax = gca;  % get the current axis
ax.Clipping = 'off';  % turn of clipping when magnifying
%% Lighting for curved membrane surface
%view(10,30)
view(3);
%shading interp
material shiny 
lightangle(30,-30)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.4;
h.DiffuseStrength = 0.7;
h.SpecularStrength = 0.7;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
camproj('perspective')
camlight(120,30)
camlight(90,60)
camlight(180,-60)
end