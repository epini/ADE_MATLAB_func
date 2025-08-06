%% Test Anisotropic Diffusive Equation MATLAB functions

% This MATLAB script runs a test of the Anisotropic Diffusive Equation
% functions used to evaluate the transmittance/reflectance for an 
% anisotropic slab of thickness L [μm].
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform, mua [1/μm], or absent for the
% steady-state intensity profiles.
% t is an array of times [ps], x and y are an array of positions [μm],
% while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity
% gaussian distribution at t = 0 along x and y, defined only for Rxyt and Txyt.

% Note that when discretizing a continuous function the binning width in
% space and time must be accounted for (e.g. mean(diff(t)).

%% set parameters
% All parameters are set in μm and ps, but you can assume they are mm and ns
% if preferred due to the rescaling property of the diffusive equation.

% if computation is too slow, consider reducing time and space bins or the
% number of virtual sources considered in the expansion M in the functions

L = 1000;
n_in = 1.3;
n_ext = 1;
t = 0:500;
mua = 3e-5;
lx = 30;
ly = 10;
lz = 20;

%% Total Time-Resolved reflectance/transmittance
% this corresponds to integrating all of the reflectance/transmittance
% spatially over the whole slab

Rt = Rt_ADE(t, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(t));
Tt = Tt_ADE(t, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(t));

figure(1), hold on, grid on, box on
plot(t, Rt, 'r')
plot(t, Tt, 'b')
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 16)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 16)
legend('Total Reflectance', 'Total Transmittance','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 1e-8 1])

%% Time-Resolved reflectance/transmittance at different locations

sx = 20; % set standard deviation a t = 0 along x [μm]
sy = 20; % set standard deviation a t = 0 along y [μm]

x = linspace(400, 2000, 5); % define positions of collection [μm]
y = linspace(400, 2000, 5); % define positions of collection [μm]

figure(2), hold on
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
    Rxyt = squeeze(Rxyt_ADE(x(i), 0, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)) * mean(diff(t));
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
    Rxyt = squeeze(Rxyt_ADE(0, y(i), t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)) * mean(diff(t));
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
    Txyt = squeeze(Txyt_ADE(x(i), 0, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)) * mean(diff(t));
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
    Txyt = squeeze(Txyt_ADE(0, y(i), t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)) * mean(diff(t));
    plot(t, Txyt)
    legend_entries{i} = sprintf('(0, %g)', y(i));
end
set(gca, 'yscale', 'log')
ylabel('Intensity [a.u.]','interpreter','latex', 'Fontsize', 14)
xlabel('$t$ [ps]','interpreter','latex', 'Fontsize', 14)
axis([0 max(t) 5e-16 1e-8])
legend(legend_entries, 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 12)


%% Time- and Space-Resolved reflectance
% transmittance is identical after the ballistic transient when the frames are normalized

x = -500:10:500; % define grid for image
y = -500:10:500;

t2 = linspace(0.5, 15.5, 16); % define shorter times, adjust the number of tiles if the number of times is changed

Rxyt = Rxyt_ADE(x, y, t2, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)*mean(diff(t2))*mean(diff(x))*mean(diff(y));
Txyt = Txyt_ADE(x, y, t2, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)*mean(diff(t2))*mean(diff(x))*mean(diff(y));

figure(3) % these frames are in linear scale and normalized to their maximum
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(4,4);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

for i = 1:length(t2)
    figure(3), nexttile
    imagesc(Rxyt(:,:,i));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    axis equal tight
    title(strcat(num2str(t2(i)), ' ps'),'interpreter','latex', 'Fontsize', 16)
end

%% Steady-state reflectance/transmittance

Rxy = Rxy_ADE(x, y, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));
Txy = Txy_ADE(x, y, L, n_in, n_ext, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));

figure(4), hold on % these frames are in logarithmic scale, mind the different colorscales for Rxy and Txy
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
