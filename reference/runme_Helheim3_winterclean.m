clear all;close all

% Script to run SHAKTI stand-alone model on Helheim Glacier
addpath Code/

steps = [1:4];
clustername = oshostname();

org=organizer('repository',['./Models'],'prefix',['Model_Helheim_'],'steps',steps); clear steps;

if perform(org,'Mesh'),% {{{
    
    % Initialize mesh as uniform with typical edge length of 500 m
    md=triangle(model,['Exp/Helheim14.exp'],500);
    
    % Load surface velocity
    [velx vely]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
    vel  = sqrt(velx.^2+vely.^2);
    
    % Refine mesh using surface velocities as metric
    md=bamg(md,'hmin',200,'hmax',10000,'field',vel,'err',5);
    [md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
    md.mesh.epsg=3413;
    
    savemodel(org,md);
end %}}}

if perform(org,'Param'),% {{{
    
    md=loadmodel(org,'Mesh');
    md=setflowequation(md,'SSA','all');
    
    md=setmask(md,'','');
    md=parameterize(md,'Par/Greenland.par');
    md.miscellaneous.name = 'Helheim';
    
    md.initialization.pressure = 0.5*md.materials.rho_ice*md.constants.g*md.geometry.thickness;
    md.initialization.temperature = 273*ones(md.mesh.numberofvertices,1);
    
    % Set sliding velocity as portion of velocity and smooth
    ub_frac=1.0;
    md.initialization.vx=ub_frac*md.initialization.vx;
    md.initialization.vy=ub_frac*md.initialization.vy;
    md.initialization.vx=averaging(md,md.initialization.vx,10);
    md.initialization.vy=averaging(md,md.initialization.vy,10);
    md.initialization.vel=ub_frac*md.initialization.vel;
    md.initialization.vel=averaging(md,md.initialization.vel,10);
    
    %flow law initialization
    disp('      creating flow law parameters (assume ice is at 0C)');
    md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
    md.materials.rheology_B = cuffey(273.15-0)*ones(md.mesh.numberofvertices,1); % Use flow law parameter for 0 deg ice
    
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
    md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
    md.basalforcings.geothermalflux = 0.05*ones(md.mesh.numberofvertices,1);
    
    % HYDROLOGY SPECIFIC PARAMETERIZATION:
    % Change hydrology class to SHAKTI model
    md.hydrology=hydrologyshakti();
    
    % Define distributed englacial input to the subglacial system (m/yr)
    md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);
    
    % Define initial water head such that water pressure is 50% of ice overburden pressure
    md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;
    
    % Initial subglacial gap height of 0.01m everywhere
    md.hydrology.gap_height = 0.01*ones(md.mesh.numberofelements,1);
    
    % Typical bed bump bump spacing
    md.hydrology.bump_spacing = 1.0*ones(md.mesh.numberofelements,1);
    
    % Typical bed bump height
    md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);
    
    % Initial Reynolds number (start at Re=1000 everywhere)
    md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);
    
    % Glacier front b.c.
    md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
    terminus=ContourToNodes(md.mesh.x,md.mesh.y,'Exp/terminus13.exp',1);
    pos=find(terminus & md.mesh.vertexonboundary);
    md.hydrology.spchead(pos)=0;
    
    md.hydrology.moulin_input = zeros(md.mesh.numberofvertices,1); % No moulin inputs
    md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1); % No-flux b.c. on boundary except outflow
    
    savemodel(org,md);
end%}}}


if perform(org,'Hydro'),% {{{
    disp('	Step 4: Hydro');
    md=loadmodel(org,'Param');
    
    % Use frictionshakti for stand-alone SHAKTI
    md.friction=frictionshakti;
    
    % Define friction coefficient inverted from end of winter spin-up
    load Data/friction_coefficient_Nfinal.mat
    md.friction.coefficient=friction_coefficient;
    fast=find(md.initialization.vel>=6000);
    md.friction.coefficient(fast)=min(md.friction.coefficient(fast));
    
    md.transient=deactivateall(md.transient);
    md.transient.ishydrology=1;
    md.inversion.iscontrol=0;
    
    % Request basal meltrate and components as outputs
    md.transient.requested_outputs={'Dummy','EsaEmotion','EsaNmotion','EsaUmotion'};
    
    % Specify that you want to run the model on your current computer
    % Change the number of processors according to your machine (here np=4)
    md.cluster=generic('np',4);
    
    % Define the time stepping scheme:
    md.timestepping.time_step=1800/md.constants.yts; % Time step (in years)
    md.timestepping.final_time=1/365;
    md.settings.output_frequency=1;
    
    md.stressbalance.restol=0.05;
    md.stressbalance.reltol=0.05;
    md.stressbalance.abstol=NaN;
    md.settings.solver_residue_threshold=NaN;
    
    md.verbose.solution=1;
    md=solve(md,'Transient');
    
    f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
    Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number
    
    savemodel(org,md);
end
