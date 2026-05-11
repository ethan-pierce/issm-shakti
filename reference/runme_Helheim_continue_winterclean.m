clear all;close all
load Models/Helheim_winter_A0_1d.mat

% Starting conditions
md.hydrology.head=md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end).HydrologyBasalFlux./1.787e-6;

md.hydrology.englacial_input(:)=0.0;

% Time-stepping
md.timestepping.time_step=1800/md.constants.yts; % Time step (in years)
md.timestepping.final_time=365/365;
md.settings.output_frequency=48;

md.transient.requested_outputs={'Dummy','EsaEmotion','EsaNmotion','EsaUmotion'};

md.cluster=generic('np',20);
md.cluster.interactive=0;
md.settings.waitonlock=0;
md.verbose.solution=1;
md=solve(md,'Transient');

f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
f(f<0)=0;
Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

%md=loadresultsfromcluster(md,'runtimename','Helheim-08-14-2022-14-06-28-16925');
%save('Models/Helheim_winter_A0_1yr','md','f','Re')

