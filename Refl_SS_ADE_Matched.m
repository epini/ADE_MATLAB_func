function Rxy = Refl_SS_ADE_Matched(x, y, L, n_matched, lx, ly, lz)

% this function returns the space resolved steady-state reflectance R(x,y)
% for a non-absorbing anisotropic slab of thickness L [um].
% The refractive index is matched with the environment.
% lx, ly and lz are scalars in microns.

v=299.7924589/n_matched;

if lx == lz && lx == ly
    Dx = lx*v/3;
    Dy = Dx;
    Dz = lz*v/3;
    ze = 2*lx/3;
else
    mux=1/lx;
    muy=1/ly;
    muz=1/lz;

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
    
    Cfun = @(chi,phi) (chi.^2)./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2)).^2;
    Bfun = @(chi,phi)  chi./(mux.*(1 - chi.^2).*((cos(phi)).^2) + muy.*(1 - chi.^2).*((sin(phi)).^2) + muz.*(chi.^2));
    C = integral2(Cfun, 0, 1, 0, 2*pi);
    B = integral2(Bfun, 0, 1, 0, 2*pi);
    ze=C/B;
end

D = (Dx*Dy*Dz)^(1/3);
z0 = lz;

Rxy = zeros(length(x),length(y));

M = 5000; %number of iterations
for m = -M:M
    z3 = -2*m*L - 4*m*ze - z0;
    z4 = -2*m*L - (4*m - 2)*ze + z0;
    Rxy = Rxy + z3.*(1/Dz*z3^2 + 1/Dx.*(x.^2).' + 1/Dy.*y.^2).^(-3/2) - z4.*(1/Dz*z4^2 + 1/Dx.*(x.^2).' + 1/Dy.*y.^2).^(-3/2);
end

Rxy = -D^(-3/2)*Rxy/4/pi;

end