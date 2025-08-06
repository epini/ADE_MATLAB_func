%% Example of space-resolved reflectance/transmittance at different positions

% This MATLAB script runs an example of the Anisotropic Diffusive Equation
% functions used to evaluate the Space-Resolved transmittance/reflectance 
% for an anisotropic slab of thickness L [μm].
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform and isotropic, mua [1/μm].
% x and y are an array of positions [μm],
% while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity
% gaussian distribution at t = 0 along x and y.

% Note that when discretizing a continuous function the binning width in
% space must be accounted for (e.g. mean(diff(x)).

% Author:       Ernesto Pini
% Affiliation:  Department of Physics and Astronomy, Università di Firenze
% Email:        pinie@lens.unifi.it

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
x = -500:10:500; % define grid for frames
y = -500:10:500;
sx = 10; % set standard deviation a t = 0 along x [μm]
sy = 10; % set standard deviation a t = 0 along y [μm]

%% run test and plot results

Rxy = Rxy_ADE(x, y, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));
Txy = Txy_ADE(x, y, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));

figure(1), hold on % these frames are in logarithmic scale, mind the different colorscales for Rxy and Txy
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(2,2);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

nexttile
imagesc(x, y, Rxy.')
axis equal tight
colorbar
set(gca,'ColorScale','log')
axis([min(x) max(x) min(y) max(y)])
ylabel('$y$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
xlabel('$x$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
title('Steady-state Reflectance','interpreter','latex', 'Fontsize', 14)

nexttile
imagesc(x, y, Txy.')
axis equal tight
colorbar
set(gca,'ColorScale','log')
axis([min(x) max(x) min(y) max(y)])
ylabel('$y$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
xlabel('$x$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
title('Steady-state Transmittance','interpreter','latex', 'Fontsize', 14)

nexttile, grid on, box on
contour(x, y, Rxy.', logspace(0,-8,33), 'k')
axis equal tight
axis([min(x) max(x) min(y) max(y)])
ylabel('$y$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
xlabel('$x$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
title('Steady-state Reflectance','interpreter','latex', 'Fontsize', 14)
legend('$10^{0.25\times n}$','interpreter','latex', 'Fontsize', 12)

nexttile, grid on, box on
contour(x, y, Txy.', logspace(0,-8,33), 'k')
axis equal tight
axis([min(x) max(x) min(y) max(y)])
ylabel('$y$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
xlabel('$x$ [$\mu$m]','interpreter','latex', 'Fontsize', 14)
title('Steady-state Transmittance','interpreter','latex', 'Fontsize', 14)
legend('$10^{0.25\times n}$','interpreter','latex', 'Fontsize', 12)
