%% Example of total reflectance/transmittance

% This MATLAB script runs an example of the Anisotropic Diffusive Equation
% functions used to evaluate the total transmittance/reflectance
% for an anisotropic slab of thickness L [μm].
% This corresponds to integrating all of the reflectance/transmittance
% spatially and temporally over the whole slab surface.
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
% Absorption is considered to be uniform and isotropic, mua [1/μm].
% lx, ly and lz are scalars [μm].

% Author:       Ernesto Pini
% Affiliation:  Department of Physics and Astronomy, Università di Firenze
% Email:        pinie@lens.unifi.it

%% set parameters
% All parameters are set in μm, but you can assume they are mm
% if preferred due to the rescaling property of the diffusive equation.

% if computation is too slow, consider reducing the number of virtual 
% sources considered in the expansion M in the functions

L = 100;
n_in = 1.3;
n_ext = 1;
mua = 3e-5;
lx = 30;
ly = 10;
lz = 20;

%% run test and plot results

R = R_ADE(L, n_in, n_ext, lx, ly, lz, mua);
T = T_ADE(L, n_in, n_ext, lx, ly, lz, mua);

disp('Reflectance')
disp(R)
disp('Transmittance')
disp(T)
disp('Absorbance')
disp(1-T-R)