%% Test Anisotropic Diffusive Equation MATLAB functions

% This MATLAB script runs a test of the Anisotropic Diffusive Equation
% functions used to evaluate the transmittance/reflectance for an 
% anisotropic slab of thickness L [μm].
% The refractive index is matched with the environment.
% Absorption is considered to be uniform, mua [1/μm], or absent for the
% steady-state intensity profiles.
% t is an array of times [ps], x and y are an array of positions [μm],
% while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity
% gaussian distribution at t = 0 along x and y, defined only for Rxyt and Txyt.

% Note that when discretizing a continuous function the binning width in
% space and time must be accounted for (e.g. mean(diff(t)).

%% set parameters
% if computation is too slow, consider reducing time and space bins or the
% number of virtual sources considered in the expansion M in the functions

L = 1000;
n_matched = 1.3;
t = 0:500;
x = -500:10:500;
y = -500:10:500;
mua = 3e-5;
lx = 30;
ly = 10;
lz = 20;

%% Total Time-Resolved reflectance/transmittance

Rt = Rt_ADE_matched(t, L, n_matched, lx, ly, lz, mua)*mean(diff(t));
Tt = Tt_ADE_matched(t, L, n_matched, lx, ly, lz, mua)*mean(diff(t));

figure(1), hold on, grid on, box on
plot(t, Rt, 'r')
plot(t, Tt, 'b')
set(gca, 'yscale', 'log')
ylabel('$\displaystyle$ Intensity (a.u.)','interpreter','latex', 'Fontsize', 16)
xlabel('$\displaystyle t$ (ps)','interpreter','latex', 'Fontsize', 16)
legend('Total Reflectance', 'Total Transmittance','interpreter','latex', 'Fontsize', 14)
axis([0 length(t)*mean(diff(t)) 1e-8 1])

%% Time- and Space-Resolved reflectance
% transmittance is identical after the ballistic transient when the frames are normalized

t2 = linspace(0.5, 15.5, 16); % define shorter times, adjust the number of tiles if the number of times is changed

sx = 20; % set standard deviation a t = 0 along x [μm]
sy = 20; % set standard deviation a t = 0 along x [μm]

Rxyt = Rxyt_ADE_matched(x, y, t2, L, n_matched, lx, ly, lz, sx, sy, mua)*mean(diff(t2))*mean(diff(x))*mean(diff(y));
Txyt = Txyt_ADE_matched(x, y, t2, L, n_matched, lx, ly, lz, sx, sy, mua)*mean(diff(t2))*mean(diff(x))*mean(diff(y));

figure(2) % these frames are in linear scale
x0 = 500;
y0 = 200;
width = 500;
height = 500;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(4,4);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

for i = 1:length(t2)
    figure(2), nexttile
    imagesc(Rxyt(:,:,i));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    axis equal tight
    title(strcat('$\displaystyle$', num2str(t2(i)), ' ps'),'interpreter','latex', 'Fontsize', 10)
end

%% Steady-state reflectance/transmittance

Rxy = Rxy_ADE_matched(x, y, L, n_matched, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));
Txy = Txy_ADE_matched(x, y, L, n_matched, lx, ly, lz, mua)*mean(diff(x))*mean(diff(y));

figure(3) % these frames are in logarithmic scale, mind the different colorbar for Rxy and Txy
x0 = 350;
y0 = 60;
width = 800;
height = 800;
set(gcf,'position',[x0,y0,width,height])
tile = tiledlayout(2,2);
tile.Padding = 'compact';
tile.TileSpacing = 'compact';
colormap(parula(256))

figure(3), nexttile
imagesc(x, y, Rxy)
axis equal tight
colorbar
set(gca,'ColorScale','log')
axis([min(x) max(x) min(y) max(y)])
ylabel('$\displaystyle y$ (micron)','interpreter','latex', 'Fontsize', 14)
xlabel('$\displaystyle x$ (micron)','interpreter','latex', 'Fontsize', 14)
title('$\displaystyle$ Steady-state Reflectance','interpreter','latex', 'Fontsize', 14)

figure(3), nexttile
imagesc(x, y, Txy)
axis equal tight
colorbar
set(gca,'ColorScale','log')
axis([min(x) max(x) min(y) max(y)])
ylabel('$\displaystyle y$ (micron)','interpreter','latex', 'Fontsize', 14)
xlabel('$\displaystyle x$ (micron)','interpreter','latex', 'Fontsize', 14)
title('$\displaystyle$ Steady-state Transmittance','interpreter','latex', 'Fontsize', 14)

figure(3), nexttile, grid on, box on
contour(x, y, Rxy, logspace(0,-8,33), 'k')
axis equal tight
axis([min(x) max(x) min(y) max(y)])
ylabel('$\displaystyle y$ (micron)','interpreter','latex', 'Fontsize', 14)
xlabel('$\displaystyle x$ (micron)','interpreter','latex', 'Fontsize', 14)
title('$\displaystyle$ Steady-state Reflectance','interpreter','latex', 'Fontsize', 14)
legend('$\displaystyle$ $10^{0.25\times n}$','interpreter','latex', 'Fontsize', 12)

figure(3), nexttile, grid on, box on
contour(x, y, Txy, logspace(0,-8,33), 'k')
axis equal tight
axis([min(x) max(x) min(y) max(y)])
ylabel('$\displaystyle y$ (micron)','interpreter','latex', 'Fontsize', 14)
xlabel('$\displaystyle x$ (micron)','interpreter','latex', 'Fontsize', 14)
title('$\displaystyle$ Steady-state Transmittance','interpreter','latex', 'Fontsize', 14)
legend('$\displaystyle$ $10^{0.25\times n}$','interpreter','latex', 'Fontsize', 12)
