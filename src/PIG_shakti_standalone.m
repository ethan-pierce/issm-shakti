% Load ISSM tools
addpath("/Users/f007s79/Documents/repos/ISSM-binaries/bin");
addpath("/Users/f007s79/Documents/repos/ISSM-binaries/lib");
issmversion;

%%%%%%%%%%%%%%%%
% Load the model
%%%%%%%%%%%%%%%%
load ../models/PIG_SHAKTI_spinup md;

md.hydrology.head = md.results.TransientSolution(end).HydrologyHead;
md.hydrology.gap_height = md.results.TransientSolution(end).HydrologyGapHeight;
md.hydrology.reynolds = md.results.TransientSolution(end).HydrologyBasalFlux ./ 1.787e-6;
md.hydrology.reynolds(md.hydrology.reynolds < 100) = 100;

md.hydrology.englacial_input(:) = 0.0;
md.hydrology.moulin_input(:) = 0.0;
md.hydrology.neumannflux(:) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%
% Set up new simulation
%%%%%%%%%%%%%%%%%%%%%%%
md.timestepping.time_step = 1000 / md.constants.yts;

md.timestepping.final_time = 180 / 365;
md.settings.output_frequency = 40;

md.cluster = generic('np', 10);
md.verbose.solution = 1;
md = solve(md, 'Transient');

Pi = md.materials.rho_ice * md.constants.g * md.geometry.thickness;
Pw = md.materials.rho_freshwater * md.constants.g * (md.results.TransientSolution(end).HydrologyHead - md.geometry.base);
N = Pi - Pw;

% Plot the results
plotmodel(md, 'data', N);
exportgraphics(gcf, '../figures/standalone/PIG_shakti_N.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.results.TransientSolution(end).HydrologyHead);
exportgraphics(gcf, '../figures/standalone/PIG_shakti_head.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.results.TransientSolution(end).HydrologyGapHeight);
exportgraphics(gcf, '../figures/standalone/PIG_shakti_gap_height.png', 'Resolution', 300);
close all;

plotmodel(md, 'data', md.results.TransientSolution(end).HydrologyBasalFlux);
exportgraphics(gcf, '../figures/standalone/PIG_shakti_basal_flux.png', 'Resolution', 300);
close all;

save ../models/PIG_SHAKTI_standalone md;
