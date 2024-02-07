function Rxyt = Refl_TSR_ADE_Matched(x, y, t, L, n_matched, lx, ly, lz, sx, sy, mua)

% this function returns the time- and space-resolved Reflectance R(x,y,t)
% for an anisotropic slab of thickness L [um].
% The refractive index is matched with the environment. Absorption is considered to be uniform, mua [1/um].
% t is an array of times in ps, while lx, ly and lz are scalars in microns.
% sx and sy are the initial width in microns of the intensity gaussian distribution at t = 0 along x and y.

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

R = zeros(size(t));
Rxyt = zeros(length(x),length(y),length(t));

M = 10000; %number of iterations
for m = -M:M
    z3 = -2*m*L - 4*m*ze - z0;
    z4 = -2*m*L - (4*m - 2)*ze + z0;
    R = R + (z3*exp(-(z3)^2./(4*Dz*t))-z4*exp(-(z4)^2./(4*Dz*t)));
end

j=1;
for i=t
    Rxyt(:,:,j) = -exp(-x.^2./(sx^2+4*Dx*i)).*(exp(-y.^2./(sy^2+4*Dy*i))).'./((2*(4*pi)^(3/2)*(4*i)^(5/2)*(Dx*Dy*Dz)^(1/2))).*R(j).'.*exp(-v*i*mua)/2;
    j=j+1;
end
end