%% Example of total time-resolved reflectance/transmittance

% This MATLAB script runs an example of the Anisotropic Diffusive Equation
% functions used to evaluate the total time-resolved transmittance/reflectance
% for an anisotropic slab of thickness L [μm].
% This corresponds to integrating all of the time-resolved
% reflectance/transmittance spatially over the whole slab surface.
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform and isotropic, mua [1/μm].
% t is an array of times [ps], while lx, ly and lz are scalars [μm].

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

%% run test and plot results

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