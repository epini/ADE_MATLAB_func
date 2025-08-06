function [Dx, Dy, Dz] = D_Tensor_ADE(n_in, lx, ly, lz)
% D_TENSOR_ADE Diffusive rate tensor from scattering mean free paths
%
% Brief: This function returns the (diagonal) diffusive rate tensor
% elements Dx, Dy, Dz starting from the scattering mean free paths
% lx, ly and lz. Quantities can be expressed either in μm and ps or in
% mm and ns.
%
% Inputs:
%    lx - scattering mean free path along x [μm] or [mm]
%    ly - scattering mean free path along y [μm] or [mm]
%    lz - scattering mean free path along z [μm] or [mm]
%    n_in - refractive index of the diffusive medium
%
% Outputs:
%    Dx - diffusive rate along x [μm^2/ps] or [mm^2/ns]
%    Dy - diffusive rate along y [μm^2/ps] or [mm^2/ns]
%    Dz - diffusive rate along z [μm^2/ps] or [mm^2/ns]
%
% See also: Test_function.m

% Author:       Ernesto Pini
% Affiliation:  Department of Physics and Astronomy, Università di Firenze
% Email:        pinie@lens.unifi.it

v = 299.7924589/n_in;

if lx == lz && lx == ly % the isotropic case is treated separately

    Dx = lx*v/3;
    Dy = Dx;
    Dz = Dx;

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
    Dy = v/(4*pi*lavg)*integral2(Dyfun, -1, 1, 0, 2*pi);
    Dz = v/(4*pi*lavg)*integral2(Dzfun, -1, 1, 0, 2*pi);

end

end