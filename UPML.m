
Nx = 80;
Ny = 80;
Npml = [10,10]
sigx = zeros(2*Nx, 2*Ny);
dt = 0.1;
for nx = 1:2*Npml(1)
   nx1 = 2*Npml(1)-nx+1;
   sigx(nx1,:) = 0.5*dt*(nx/2/Npml(1))^3;
end

