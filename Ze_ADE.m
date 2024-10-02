function ze = Ze_ADE(n_in, n_ext, lx, ly, lz)
% ZE_ADE Extrapolated length for a slab geometry
%
% Brief: This function returns the extrapolated length from an anisotropic
% slab.
% xy is the slab plane, while z is the direction of incidence of the pencil beam.
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered.
%
% Inputs:
%    n_in - refractive index of the diffusive medium
%    n_ext - refractive index of the external medium
%    lx - scattering mean free path along x [μm]
%    ly - scattering mean free path along y [μm]
%    lz - scattering mean free path along z [μm]
%
% Outputs:
%    ze - extrapolated length
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

    if n == 1
        ze = 2*lx/3;
    else
        A = (1 + 3*integral(Cfun, 0, 1))./(1 - 2*integral(Bfun, 0, 1));
        ze = 2*A*lx/3;
    end

else % anisotropic case
    mux = 1/lx;
    muy = 1/ly;
    muz = 1/lz;

    Rfun = @(chi) (1/2).*(abs((n.*chi - sqrt(1 - (1 - chi.^2).*n^2))./(n*chi + sqrt(1 - (1 - chi.^2).*n^2))).^2 + abs((chi - n.*sqrt(1 - (1 - chi.^2).*n^2))./(chi + n.*sqrt(1 - (1 - chi.^2).*n^2))).^2);
    Cfun = @(chi,phi) (chi.^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;
    Bfun = @(chi,phi)  chi./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2));
    Xfun = @(chi,phi) Bfun(chi,phi).*Rfun(chi);
    Yfun = @(chi,phi) Cfun(chi,phi).*Rfun(chi);

    C = integral2(Cfun, 0, 1, 0, 2*pi);
    B = integral2(Bfun, 0, 1, 0, 2*pi);

    if n == 1
        ze = C/B;
    else
        X = integral2(Xfun, 0, 1, 0, 2*pi);
        Y = integral2(Yfun, 0, 1, 0, 2*pi);
        ze = (C + Y)/(B - X);
    end
end

end