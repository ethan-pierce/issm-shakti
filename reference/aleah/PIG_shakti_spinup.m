%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize and refine the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
domain = ['PIGOutline.exp'];
hinit = 5000; %10000
hmax = 5000; %10000
hmin = 500; %2500
gradation = 1.7;
err = 8;
md = bamg(model, 'domain', domain, 'hmax', hinit);

velocity = 'Antarctica_ice_velocity.nc';
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

%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterize the model
%%%%%%%%%%%%%%%%%%%%%%%%
md = setflowequation(md, 'SSA', 'all');
md = setmask(md, '', '');
paramfile = 'PIG.par';
md = parameterize(md, paramfile);
md.miscellaneous.name = 'PIG';

md.initialization.pressure = 0.5 * md.materials.rho_ice * md.constants.g * md.geometry.thickness;
md.initialization.temperature = (273-5) * ones(md.mesh.numberofvertices, 1); % Ethan 273
ub_frac = 1.0;
md.initialization.vx = ub_frac * md.initialization.vx;
md.initialization.vy = ub_frac * md.initialization.vy;
md.initialization.vx = averaging(md, md.initialization.vx, 10);
md.initialization.vy = averaging(md, md.initialization.vy, 10);
md.initialization.vel = ub_frac * md.initialization.vel;
md.initialization.vel = averaging(md, md.initialization.vel, 10);

md.materials.rheology_n = 3 * ones(md.mesh.numberofelements, 1);
md.materials.rheology_B = cuffey(273.15 - 5) * ones(md.mesh.numberofvertices, 1); % Ethan temperate

md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices, 1);
md.basalforcings.geothermalflux = 0.05 * ones(md.mesh.numberofvertices, 1);

% plotmodel(md, 'data', md.geometry.thickness);
% exportgraphics(gcf, '../figures/initial/PIG_thickness.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.geometry.surface);
% exportgraphics(gcf, '../figures/initial/PIG_surface.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.geometry.base);
% exportgraphics(gcf, '../figures/initial/PIG_base.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.mask.ocean_levelset);
% exportgraphics(gcf, '../figures/initial/PIG_ocean_levelset.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.mask.ice_levelset);
% exportgraphics(gcf, '../figures/initial/PIG_ice_levelset.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.initialization.vel);
% exportgraphics(gcf, '../figures/initial/PIG_initial_vel.png', 'Resolution', 300);
% close all;

%%%%%%%%%%%%%%%%%%%%%
% Parameterize SHAKTI
%%%%%%%%%%%%%%%%%%%%%
md.hydrology = hydrologyshakti();
md.hydrology.englacial_input = 0.0 * ones(md.mesh.numberofvertices, 1);
md.hydrology.head = 0.5 * md.materials.rho_ice / md.materials.rho_freshwater * md.geometry.thickness + md.geometry.base;
md.hydrology.gap_height = 0.01 * ones(md.mesh.numberofelements, 1);
md.hydrology.bump_spacing = 1.0 * ones(md.mesh.numberofelements, 1);
md.hydrology.bump_height = 0.0 * ones(md.mesh.numberofelements, 1);
md.hydrology.reynolds = 1000 * ones(md.mesh.numberofelements, 1);

% Glacier front b.c.
md.hydrology.spchead = NaN(md.mesh.numberofvertices, 1);
% Option 1)
ocean = find(md.mask.ocean_levelset <= 0);
md.hydrology.spchead(ocean) = 0;

% Option 2)
% md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);
% pos = find(md.mask.ice_levelset >= -50); % Find everywhere within 50m inland
% md.hydrology.spchead(pos) = 0;

% Option 3)
% md.mask.ocean_levelset(:) = 1;
% md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);
% pos = find(md.mask.ice_levelset >= -1000 | md.mask.ocean_levelset < 0);
% md.hydrology.spchead(pos) = 0;

md.hydrology.moulin_input = zeros(md.mesh.numberofvertices, 1); 
md.hydrology.neumannflux = zeros(md.mesh.numberofelements, 1);

% plotmodel(md, 'data', md.hydrology.spchead);
% exportgraphics(gcf, '../figures/initial/PIG_shakti_spchead.png', 'Resolution', 300);
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SHAKTI in decoupled mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% md.friction = frictionshakti;
% load PIG_inversion_friction.mat
% md.friction.coefficient = fcoeff;

% ANS
% md.friction=();
md.friction.coefficient(:)=0;

md.transient = deactivateall(md.transient);
md.transient.ishydrology = 1;
md.inversion.iscontrol = 0;
md.cluster = generic('np', 4);
md.timestepping.time_step = 3600 / md.constants.yts; % Ethan 1800
md.timestepping.final_time = 10 / 365;
md.settings.output_frequency = 1;
md.stressbalance.restol = 0.05;
md.stressbalance.reltol = 0.05;
md.stressbalance.abstol = NaN;
md.settings.solver_residue_threshold = NaN;
md.verbose.solution = 1;
md = solve(md, 'Transient');

f = md.materials.rho_freshwater ./ (md.materials.rho_ice .* md.geometry.thickness) ...
    .* (md.results.TransientSolution(end).HydrologyHead - md.geometry.base); 
Re = abs(md.results.TransientSolution(end).HydrologyBasalFlux) ./ 1.787e-6;
Pi = md.materials.rho_ice * md.constants.g * md.geometry.thickness;
Pw = md.materials.rho_freshwater * md.constants.g * (md.results.TransientSolution(end).HydrologyHead - md.geometry.base);
N = Pi - Pw;

% plotmodel(md, 'data', Pw);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_Pw.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.results.TransientSolution(end).HydrologyBasalFlux);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_basal_flux.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', md.results.TransientSolution(end).HydrologyGapHeight);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_gap_height.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', N);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_N.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', f);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_fraction.png', 'Resolution', 300);
% close all;

% plotmodel(md, 'data', Re);
% exportgraphics(gcf, '../figures/spinup/PIG_shakti_reynolds.png', 'Resolution', 300);
% close all;

% save ../models/PIG_SHAKTI_spinup md;