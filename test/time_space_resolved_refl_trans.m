%% Test Anisotropic Diffusive Equation MATLAB functions for
%% Time- and Space-Resolved reflectance/transmittance at different positions

% This MATLAB script runs a test of the Anisotropic Diffusive Equation
% functions used to evaluate the Time- and Space-Resolved
% transmittance/reflectance for an anisotropic slab of thickness L [μm].
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform, mua [1/μm], or absent for the
% steady-state intensity profiles.
% t is an array of times [ps], x and y are an array of positions [μm],
% while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity
% gaussian distribution at t = 0 along x and y, defined only for Rxyt and Txyt.

% Note that when discretizing a continuous function the binning width in
% time and space must be accounted for (e.g. mean(diff(t)).

%% set parameters
% All parameters are set in μm and ps, but you can assume they are mm and ns
% if preferred due to the rescaling property of the diffusive equation.

% if computation is too slow, consider reducing time and space bins or the
% number of virtual sources considered in the expansion M in the functions

L = 1000;
n_in = 1.3;
n_ext = 1;
mua = 3e-5;
lx = 30;
ly = 10;
lz = 20;
t = linspace(0.2, 12.2, 16);
x = -500:10:500; % define grid for frames
y = -500:10:500;
sx = 10; % set standard deviation a t = 0 along x [μm]
sy = 10; % set standard deviation a t = 0 along y [μm]

%% run test and plot results
% reflectance and transmittance frames are identical after the ballistic
% transient if the frames are normalized, only the amplitude is different

Rxyt = Rxyt_ADE(x, y, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)*mean(diff(t))*mean(diff(x))*mean(diff(y));
Txyt = Txyt_ADE(x, y, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)*mean(diff(t))*mean(diff(x))*mean(diff(y));

figure(1) % reflectance frames
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(4,4);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

for i = 1:length(t)
    figure(1), nexttile
    imagesc(Rxyt(:,:,i));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    axis equal tight
    title(strcat(num2str(t(i)), ' ps'),'interpreter','latex', 'Fontsize', 16)
    colorbar
end

figure(2) % transmittance frames
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(4,4);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

for i = 1:length(t)
    figure(2), nexttile
    imagesc(Txyt(:,:,i));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    axis equal tight
    title(strcat(num2str(t(i)), ' ps'),'interpreter','latex', 'Fontsize', 16)
    colorbar
end