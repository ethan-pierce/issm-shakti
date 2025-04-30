% Load ISSM tools
addpath("/Users/f007s79/Documents/repos/ISSM-binaries/bin");
addpath("/Users/f007s79/Documents/repos/ISSM-binaries/lib");
issmversion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize and refine the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
domain = ['../data/PIGOutline.exp'];
hinit = 10000;
hmax = 10000;
hmin = 2500;
gradation = 1.7;
err = 8;
md = bamg(model, 'domain', domain, 'hmax', hinit);

velocity = '../data/Antarctica_ice_velocity.nc';
xmin = ncreadatt(velocity,'/','xmin');
ymax = ncreadatt(velocity,'/','ymax');
spacing	= ncreadatt(velocity,'/','spacing');
nx = double(ncreadatt(velocity,'/','nx'));
ny = double(ncreadatt(velocity,'/','ny'));
vx = double(ncread(velocity,'vx'));
vy = double(ncread(velocity,'vy'));

xmin = strtrim(xmin);
xmin = str2num(xmin(1:end-2));
ymax = strtrim(ymax);
ymax = str2num(ymax(1:end-2));
spacing = strtrim(spacing);
spacing = str2num(spacing(1:end-2));

x = xmin + (0:1:nx)'*spacing;
y = (ymax - ny * spacing) + (0:1:ny)'*spacing;

vx_obs = InterpFromGridToMesh(x, y, flipud(vx'), md.mesh.x, md.mesh.y, 0);
vy_obs = InterpFromGridToMesh(x, y, flipud(vy'), md.mesh.x, md.mesh.y, 0);
vel_obs = sqrt(vx_obs.^2 + vy_obs.^2);

md = bamg(md, 'hmax', hmax, 'hmin', hmin, 'gradation', gradation, 'field', vel_obs, 'err', err);

% plotmodel(md, 'data', 'mesh');
% exportgraphics(gcf,'../figures/PIG_mesh.png', 'Resolution', 300);
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterize the model
%%%%%%%%%%%%%%%%%%%%%%%%
md = setflowequation(md, 'SSA', 'all');
md = setmask(md, '', '');
paramfile = '../data/PIG.par';
md = parameterize(md, paramfile);
md.miscellaneous.name = 'PIG';

plotmodel(md, 'data', md.geometry.thickness);
exportgraphics(gcf, '../figures/initial/PIG_thickness.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.geometry.surface);
exportgraphics(gcf, '../figures/initial/PIG_surface.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.geometry.base);
exportgraphics(gcf, '../figures/initial/PIG_base.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.mask.ocean_levelset);
exportgraphics(gcf, '../figures/initial/PIG_ocean_levelset.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.mask.ice_levelset);
exportgraphics(gcf, '../figures/initial/PIG_ice_levelset.png', 'Resolution', 300);
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SHAKTI in decoupled mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

