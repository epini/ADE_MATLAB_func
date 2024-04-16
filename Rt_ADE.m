function Rt = Rt_ADE(t, L, n_in, n_ext, lx, ly, lz, mua)
% RT_ADE time-resolved reflectance from a turbid slab
%
% Brief: this function returns the total time-resolved reflectance R(t)
% from an anisotropic slab of thickness L [μm].
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered. Absorption is assumed to be
% uniform, mua [1/μm].
% t is an array of times [ps], while lx, ly and lz are scalars [μm].
%
% Inputs:
%    t - array of times [ps]
%    L - slab thickness [μm]
%    n_in - refractive index of the diffusive medium
%    n_ext - refractive index of the external medium
%    lx - scattering mean free path along x [μm]
%    ly - scattering mean free path along y [μm]
%    lz - scattering mean free path along z [μm]
%    mua - absorption rate [1/μm]
%
% Outputs:
%    Rt - array of total time-resolved reflectance R(t)
%
% See also: Test_function.m

% Author:       Ernesto Pini
% Affiliation:  Department of Physics and Astronomy, Università di Firenze
% Email:        pinie@lens.unifi.it

v = 299.7924589/n_in;
n = n_in/n_ext;

if lx == lz && lx == ly % the isotropic case is treat separately

    Rfun = @(chi) (1/2).*(abs((n.*chi - sqrt(1 - (1 - chi.^2).*n^2))./(n*chi + sqrt(1 - (1 - chi.^2).*n^2))).^2 + abs((chi - n.*sqrt(1 - (1 - chi.^2).*n^2))./(chi + n.*sqrt(1 - (1 - chi.^2).*n^2))).^2);
    Cfun = @(chi) (chi.^2).*Rfun(chi);
    Bfun = @(chi) chi.*Rfun(chi);
    A = (1 + 3*integral(Cfun, 0, 1))./(1 - 2*integral(Bfun, 0, 1));

    Dx = lx*v/3;
    Dy = Dx;
    Dz = Dx;
    ze = 2*A*lx/3;

else % anisotropic case
    mux = 1/lx;
    muy = 1/ly;
    muz = 1/lz;

    lavgfun = @(chi,phi) 1./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2));
    lavg = (1/4/pi)*integral2(lavgfun, -1, 1, 0, 2*pi);

    Dxfun = @(chi,phi) (1-chi.^2).*((cos(phi)).^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;
    Dyfun = @(chi,phi) (1-chi.^2).*((sin(phi)).^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;
    Dzfun = @(chi,phi) (chi.^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;

    Dx = v/(4*pi*lavg)*integral2(Dxfun, -1, 1, 0, 2*pi);
    ltx = 3*Dx/v;
    Dy = v/(4*pi*lavg)*integral2(Dyfun, -1, 1, 0, 2*pi);
    lty = 3*Dy/v;
    Dz = v/(4*pi*lavg)*integral2(Dzfun, -1, 1, 0, 2*pi);
    ltz = 3*Dz/v;

    Rfun = @(chi) (1/2).*(abs((n.*chi - sqrt(1 - (1 - chi.^2).*n^2))./(n*chi + sqrt(1 - (1 - chi.^2).*n^2))).^2 + abs((chi - n.*sqrt(1 - (1 - chi.^2).*n^2))./(chi + n.*sqrt(1 - (1 - chi.^2).*n^2))).^2);
    Cfun = @(chi,phi) (chi.^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;
    Bfun = @(chi,phi)  chi./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2));
    Xfun = @(chi,phi) Bfun(chi,phi).*Rfun(chi);
    Yfun = @(chi,phi) Cfun(chi,phi).*Rfun(chi);

    C = integral2(Cfun, 0, 1, 0, 2*pi);
    B = integral2(Bfun, 0, 1, 0, 2*pi);
    X = integral2(Xfun, 0, 1, 0, 2*pi);
    Y = integral2(Yfun, 0, 1, 0, 2*pi);
    ze = (C + Y)/(B - X);
end

D = (Dx*Dy*Dz)^(1/3);

z0 = lz;

R = zeros(size(t));

M = 10000; % number of virtual sources considered in the expansion
for m = -M:M
    z3 = - 2*m*L - 4*m*ze - z0;
    z4 = - 2*m*L - (4*m - 2)*ze + z0;
    R = R + (z3*exp(-(z3)^2./(4*Dz*t)) - z4*exp(-(z4)^2./(4*Dz*t)));
end

Rt = -(1/4)*((pi.*Dz.*(t).^3).^(-1/2)).*R.*exp(-v*t*mua);

end