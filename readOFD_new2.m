fid = fopen('ofd.out', 'r');
Title = fread(fid, 256, 'uchar');
Nx = fread(fid, 1, 'int');
Ny = fread(fid, 1, 'int');
Nz = fread(fid, 1, 'int');
Ni = fread(fid, 1, 'int');
Nj = fread(fid, 1, 'int');
Nk = fread(fid, 1, 'int');
N0 = fread(fid, 1, 'int');
NN = fread(fid, 1, 'int64');
NFreq1 = fread(fid, 1, 'int');
NFreq2 = fread(fid, 1, 'int');
NFeed = fread(fid, 1, 'int');
NPoint = fread(fid, 1, 'int');
Niter = fread(fid, 1, 'int');
Ntime = fread(fid, 1, 'int');
Solver.maxiter = fread(fid, 1, 'int');
Solver.nout = fread(fid, 1, 'int');
Dt = fread(fid, 1, 'double');
NGline = fread(fid, 1, 'int');
IPlanewave = fread(fid, 1, 'int');
Planewave.theta = fread(fid, 1, 'double');
Planewave.phi = fread(fid, 1, 'double');
Planewave.ei = fread(fid, 3, 'double');
Planewave.hi = fread(fid, 3, 'double');
Planewave.ri = fread(fid, 3, 'double');
Planewave.r0 = fread(fid, 3, 'double');
Planewave.ai = fread(fid, 1, 'double');
Planewave.pol = fread(fid, 1, 'int');

tmp = fread(fid, 1, 'int');

Xn = fread(fid, Nx + 1, 'double');
Yn = fread(fid, Ny + 1, 'double');
Zn = fread(fid, Nz + 1, 'double');

Xc = fread(fid, Nx, 'double');
Yc = fread(fid, Ny, 'double');
Zc = fread(fid, Nz, 'double');

Eiter = fread(fid, Niter, 'double');
Hiter = fread(fid, Niter, 'double');

VFeed = fread(fid, NFeed*(Solver.maxiter+1), 'double');
IFeed = fread(fid, NFeed*(Solver.maxiter+1), 'double');
VPoint = fread(fid, NPoint*(Solver.maxiter+1), 'double');
Freq1 = fread(fid,  NFreq1, 'double');
Freq2 = fread(fid,  NFreq2, 'double');
Feed = struct('dir', [], ...
    'i', [], 'j', [], 'k', [], ...
    'dx', [], 'dy', [], 'dz', [], ...
    'volt', [], 'delay', [], 'z0', []);
for i=1:NFeed
	Feed(i).dir = fread(fid,  1, 'uchar');
	Feed(i).i = fread(fid,  1, 'int');
	Feed(i).j = fread(fid,  1, 'int');
	Feed(i).k = fread(fid,  1, 'int');
	Feed(i).dx = fread(fid,  1, 'double');
	Feed(i).dy = fread(fid,  1, 'double');
	Feed(i).dz = fread(fid,  1, 'double');
	Feed(i).volt = fread(fid,  1, 'double');
	Feed(i).delay = fread(fid,  1, 'double');
	Feed(i).z0 = fread(fid,  1, 'double');
end
Zin = zeros(1, NFeed*NFreq1);
for i=1:NFeed*NFreq1
	tmpr = fread(fid, 1, 'double');
	tmpi = fread(fid, 1, 'double');
	Zin(i) = complex(tmpr, tmpi);
end
Ref = fread(fid,  NFeed*NFreq1, 'double');
Pin = zeros(2, NFeed*NFreq2);
Pin(1, :) = fread(fid, NFeed*NFreq2, 'double'); 
Pin(2, :) = fread(fid, NFeed*NFreq2, 'double'); 
Spara = zeros(1, NPoint*NFreq1);
for i=1:NPoint*NFreq1
	tmpr = fread(fid, 1, 'double');
	tmpi = fread(fid, 1, 'double');
	Spara(i) = complex(tmpr, tmpi);
end
Gline = fread(fid, NGline*2*3, 'double');
NSurface = fread(fid, 1, 'int');
%vsurface = struct('nx', [], 'ny', [], 'nz', [], 'x', [], 'y', [], 'z', [], 'ds', []);
vsurfacenx = zeros(1, NSurface);
vsurfaceny = zeros(1, NSurface);
vsurfacenz = zeros(1, NSurface);
vsurfacex = zeros(1, NSurface);
vsurfacey = zeros(1, NSurface);
vsurfacez = zeros(1, NSurface);
vsurfaceds = zeros(1, NSurface);
for i=1:NSurface
    vsurfacenx(i) = fread(fid, 1, 'double');
    vsurfaceny(i) = fread(fid, 1, 'double');
    vsurfacenz(i) = fread(fid, 1, 'double');
    vsurfacex(i) = fread(fid, 1, 'double');
    vsurfacey(i) = fread(fid, 1, 'double');
    vsurfacez(i) = fread(fid, 1, 'double');
    vsurfaceds(i) = fread(fid, 1, 'double');
end
SurfaceEx = zeros(NFreq2, NSurface);
SurfaceEy = zeros(NFreq2, NSurface);
SurfaceEz = zeros(NFreq2, NSurface);
SurfaceHx = zeros(NFreq2, NSurface);
SurfaceHy = zeros(NFreq2, NSurface);
SurfaceHz = zeros(NFreq2, NSurface);
for i=1:NFreq2
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceEx(i, j) = complex(tmpr, tmpi);
	end
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceEy(i, j) = complex(tmpr, tmpi);
	end
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceEz(i, j) = complex(tmpr, tmpi);
	end
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceHx(i, j) = complex(tmpr, tmpi);
	end
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceHy(i, j) = complex(tmpr, tmpi);
	end
	for j=1:NSurface
		tmpr = fread(fid, 1, 'double');
		tmpi = fread(fid, 1, 'double');
		SurfaceHz(i, j) = complex(tmpr, tmpi);
	end
end
cEx_r = zeros(1, NN*NFreq2);
cEx_i = zeros(1, NN*NFreq2);
cEy_r = zeros(1, NN*NFreq2);
cEy_i = zeros(1, NN*NFreq2);
cEz_r = zeros(1, NN*NFreq2);
cEz_i = zeros(1, NN*NFreq2);
cHx_r = zeros(1, NN*NFreq2);
cHx_i = zeros(1, NN*NFreq2);
cHy_r = zeros(1, NN*NFreq2);
cHy_i = zeros(1, NN*NFreq2);
cHz_r = zeros(1, NN*NFreq2);
cHz_i = zeros(1, NN*NFreq2);
for i=1:NFreq2
	cEx_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cEx_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cEy_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cEy_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cEz_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cEz_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHx_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHx_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHy_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHy_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHz_r(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
	cHz_i(NN*(i-1)+1:NN*i) = fread(fid, NN, 'float');
end
num0 = fread(fid, 1, 'uint64');

fclose(fid);
