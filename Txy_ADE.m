function Txy = Txy_ADE(x, y, L, n_in, n_ext, lx, ly, lz, mua)
% TXY_ADE space-resolved transmittance from an index-matched turbid slab
%
% Brief: this function returns the space-resolved steady-state transmittance
% T(x,y) for an anisotropic scattering slab of thickness L [μm].
% If a refractive index contrast is set, the effect of Fresnel
% reflections at the boundaries is considered. Absorption is assumed to be
% uniform, mua [1/μm].
% lx, ly and lz are scalars [μm]. The source is considered to be
% point-like.
% 
% Inputs:
%    x - array of positions [μm]
%    y - array of positions [μm]
%    L - slab thickness [μm]
%    n_in - refractive index of the diffusive medium
%    n_ext - refractive index of the external medium
%    lx - scattering mean free path along x [μm]
%    ly - scattering mean free path along y [μm]
%    lz - scattering mean free path along z [μm]
%    mua - absorption rate [1/μm]
% 
% Outputs:
%    Txy - array of space-resolved transmittance T(x,y)
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

    Dx = lx*v/3;
    Dy = Dx;
    Dz = Dx;
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
    if n == 1
        ze = C/B;
    else
        X = integral2(Xfun, 0, 1, 0, 2*pi);
        Y = integral2(Yfun, 0, 1, 0, 2*pi);
        ze = (C + Y)/(B - X);
    end
end

D = (Dx*Dy*Dz)^(1/3);
z0 = lz;

Txy = zeros(length(x),length(y));

M = 5000; % number of virtual sources considered in the expansion
for m = -M:M
    z1 = L*(1-2*m) - 4*m*ze - z0;
    z2 = L*(1-2*m) - (4*m - 2)*ze + z0;
    Txy = Txy + z1.*(z1^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(-3/2).*(1 + (mua*v)^(1/2)*(z1^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(1/2)).*exp(-((mua*v)^(1/2)*(z1^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(1/2)))...
        - z2.*(z2^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(-3/2).*(1 + (mua*v)^(1/2)*(z2^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(1/2)).*exp(-((mua*v)^(1/2)*(z2^2/Dz + (x.^2).'/Dx + y.^2/Dy).^(1/2)));
end

Txy = D^(-3/2)*Txy/4/pi;

end