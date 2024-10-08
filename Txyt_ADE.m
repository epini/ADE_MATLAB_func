function Txyt = Txyt_ADE(x, y, t, L, n_in, n_ext, lx, ly, lz, sx, sy, mua)
% TXYT_ADE time- and space-resolved transmittance from an index-matched turbid slab
%
% Brief: this function returns the time- and space-resolved transmittance 
% T(x,y,t) for an anisotropic slab of thickness L [μm].
% xy is the slab plane, while z is the direction of incidence of the pencil beam.
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered. Absorption is considered to 
% be uniform, mua [1/μm].
% t is an array of times [ps], while lx, ly and lz are scalars [μm].
% sx and sy are the initial standard deviation [μm] of the 2D intensity 
% gaussian distribution at t = 0 along x and y.
% 
% Inputs:
%    x - array of positions [μm]
%    y - array of positions [μm]
%    t - array of times [ps]
%    L - slab thickness [μm]
%    n_in - refractive index of the diffusive medium
%    n_ext - refractive index of the external medium
%    lx - scattering mean free path along x [μm]
%    ly - scattering mean free path along y [μm]
%    lz - scattering mean free path along z [μm]
%    sx - standard deviation a t = 0 along x [μm]
%    sy - standard deviation a t = 0 along y [μm]
%    mua - absorption rate [1/μm]
% 
% Outputs:
%    Txyt - 2D-array of time- and space-resolved transmittance T(x,y,t)
% 
% See also: Test_function.m

% Author:       Ernesto Pini
% Affiliation:  Department of Physics and Astronomy, Università di Firenze
% Email:        pinie@lens.unifi.it

[Dx, Dy, Dz] = D_Tensor_ADE(n_in, lx, ly, lz);
ze = Ze_ADE(n_in, n_ext, lx, ly, lz);

v = 299.7924589/n_in;
D = (Dx*Dy*Dz)^(1/3);

z0 = lz;

T = zeros(size(t));
Txyt = zeros(length(x),length(y),length(t));

M = 10000; % number of virtual sources considered in the expansion
for m = -M:M
    z1 = L*(1-2*m) - 4*m*ze - z0;
    z2 = L*(1-2*m) - (4*m - 2)*ze + z0;
    T = T + (z1*exp(-(z1)^2./(4*Dz*t)) - z2*exp(-(z2)^2./(4*Dz*t)));
end

j = 1;
for i = t
    Txyt(:,:,j) = exp(-x.^2./(2*sx^2+4*Dx*i)).*(exp(-y.^2./(2*sy^2+4*Dy*i))).'./((2*(4*pi)^(3/2)*i.^(5/2)*(Dx*Dy*Dz)^(1/2))).*T(j).'.*exp(-v*i*mua);
    j = j + 1;
end
end