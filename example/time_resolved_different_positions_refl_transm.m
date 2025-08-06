%% Example of time-Resolved reflectance/transmittance at different positions

% This MATLAB script runs an example of the Anisotropic Diffusive Equation
% functions used to evaluate the transmittance/reflectance for an 
% anisotropic slab of thickness L [μm] at different positions on the surface.
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform and isotropic, mua [1/μm].
% t is an array of times [ps], x and y are an array of positions [μm],
% while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity
% gaussian distribution at t = 0 along x and y, defined only for Rxyt and Txyt.

% Note that when discretizing a continuous function the binning width in
% time must be accounted for (e.g. mean(diff(t)).

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
t = 0:500;
x = linspace(400, 2000, 5); % define positions of collection [μm]
y = linspace(400, 2000, 5); % define positions of collection [μm]
sx = 10; % set standard deviation a t = 0 along x [μm]
sy = 10; % set standard deviation a t = 0 along y [μm]

%% run test and plot results

figure(1), hold on
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(2,2);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

% Reflectance moving along x
nexttile, hold on, grid on, box on
title('TR reflectance moving along $x$','interpreter','latex', 'Fontsize', 14)
legend_entries = cell(1, length(x));
for i = 1:length(x)
    Rxyt = squeeze(Rxyt_ADE(x(i), 0, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua))*mean(diff(t));
    plot(t, Rxyt)
    legend_entries{i} = sprintf('(%g, 0)', x(i));
end
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 14)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 5e-15 1e-8])
legend(legend_entries, 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 12)

% Reflectance moving along y
nexttile, hold on, grid on, box on
title('TR reflectance moving along $y$','interpreter','latex', 'Fontsize', 14)
legend_entries = cell(1, length(y));
for i = 1:length(y)
    Rxyt = squeeze(Rxyt_ADE(0, y(i), t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua))*mean(diff(t));
    plot(t, Rxyt)
    legend_entries{i} = sprintf('(0, %g)', y(i));
end
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 14)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 5e-16 1e-8])
legend(legend_entries, 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 12)

% Transmittance moving along x
nexttile, hold on, grid on, box on
title('TR transmittance moving along $x$','interpreter','latex', 'Fontsize', 14)
legend_entries = cell(1, length(x));
for i = 1:length(x)
    Txyt = squeeze(Txyt_ADE(x(i), 0, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua))*mean(diff(t));
    plot(t, Txyt)
    legend_entries{i} = sprintf('(%g, 0)', x(i));
end
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 14)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 5e-15 1e-8])
legend(legend_entries, 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 12)

% Transmittance moving along y
nexttile, hold on, grid on, box on
title('TR transmittance moving along $y$','interpreter','latex', 'Fontsize', 14)
legend_entries = cell(1, length(y));
for i = 1:length(y)
    Txyt = squeeze(Txyt_ADE(0, y(i), t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua))*mean(diff(t));
    plot(t, Txyt)
    legend_entries{i} = sprintf('(0, %g)', y(i));
end
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 14)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 5e-16 1e-8])
legend(legend_entries, 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 12)