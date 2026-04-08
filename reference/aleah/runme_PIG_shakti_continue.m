clear all;close all
load PIG_test4.mat

% Starting conditions
md.hydrology.head=md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end).HydrologyBasalFlux./1.787e-6;
md.hydrology.reynolds(md.hydrology.reynolds==0)=1;
md.friction.effective_pressure=md.results.TransientSolution(end).EffectivePressure;

% md.transient.isstressbalance=1; % Turn on velocity coupling
% md.initialization.vx=md.results.TransientSolution(end).Vx;
% md.initialization.vy=md.results.TransientSolution(end).Vy;
% md.initialization.vel=md.results.TransientSolution(end).Vel;

md.hydrology.englacial_input(:)=0.0;

% Time-stepping
md.timestepping.time_step=3600/md.constants.yts; % Time step (in years)
md.timestepping.final_time=20/365;
md.settings.output_frequency=1;

% md.transient.requested_outputs={'Dummy','EsaEmotion','EsaNmotion','EsaUmotion'};
md.transient.requested_outputs={'HydrologyMeltRate','HydrologyFrictionHeat','HydrologyDissipation','HydrologyPmpHeat'};

md.cluster=generic('np',8);
%md.cluster.interactive=0;
%md.settings.waitonlock=0;
md.verbose.solution=1;
md=solve(md,'Transient');

f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
f(f<0)=0;
Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

% md=loadresultsfromcluster(md,'runtimename','Helheim-07-19-2023-17-20-29-83876');
% save('Models/Helheim_SHAKTI_N_2_1yr_H10','md')

