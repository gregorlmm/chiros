%% SYSTEMLINK'S CHIRONOMID EMERGENCE MODEL
% Developed by:
% GREGORIO ALEJANDRO LÓPEZ MOREIRA MAZACOTTE (code development and implementation)
% Department of Ecohydrology and Biogeochemistry, Leibniz Institute of Freshwater Ecology and Inland Fisheries (IGB), Berlin, Germany
% gregorio.lopezmoreira@igb-berlin.de
% ALESSANDRO MANFRIN (conceptualisation)
% Institute for Environmental Sciences, RPTU - University of Kaiserslautern-Landau iES Landau, Landau, Germany
% a.manfrin@rptu.de
close ALL
clear
clc

rng('shuffle')
	
% Directory Separator
dir.dirSep = '\';

% SPECIFY DIRECTORIES
% Working directory
dir.workdir = 'C:\chirosEme_def'; % where this script is located

pth = genpath(dir.workdir);
addpath(pth)
cd(dir.workdir)

% General settings directory
dir.gensetdir = [dir.workdir '\genset'];
system(['mkdir ' dir.gensetdir]);

%% SCENARIOS: TIMING OF Bti APPLICATION AND EFFECTS OF Bti
% A) TIMING OF Bti APPLICATION
    time_app = ["I","Y","X","M"];
        % I: application is indepentend of population age structure (based on Bti application calendar)    
        % Y: application to maximise effect (when most larvae are young / young larvae peak)
        % X: application on the mix of both young and mature larvae (when total larvae population peaks)
        % M: application to minimise effect (when most larvae are mature / mature larvae peak)
        % Y, X, M -> Bti is applied every time there is a peak in the respective population curve

% B) STRENGTH OF THE EFFECT ON LARVAE MORTALITY AND TIME UNTIL EMERGENCE
	eff_size = ["m","i","s"];
        % m: mild overall effect on larvae (10% increase in mortality over 24h)
        % i: intermadiate overall effect on larvae (50% increase in mortality over 24h)
        % s: strong overall effect on larvae (90% increase in mortality over 24h)

% C) EFFECT ON TIME UNTIL EMERGENCE
	t_eme = ["A","D","NE"];   
        % A = Anticipated
        % D = Delayed
        % NE = No Effect
	        % If Bti does not affect emergence time, fGDD = 1;

% D) SEX DEPENDENCE
	sex_dep = ["SM","SF","ND"];   
        % SM = Stronger effect in Males
        % SF = Stronger effect in Females 
        % ND = No Difference
		        % If same effect on males and females
		        % fM = 1;
		        % fF = 1;
    	
%% DEFINE SCENARIOS
    % CONTROL
    sce = 1;
    scedef = table;
    scedef.ID(sce) = sce;
    scedef.time_app(sce) = "N/A";
    scedef.eff_size(sce) = "N/A";
    scedef.t_eme(sce) = "N/A";
    scedef.sex_dep(sce) = "N/A";
    
    % B.t.i. TREATMENTS
    sce = 2;
    for a = 1:length(time_app)
        for b = 1:length(eff_size)
            for c = 1:length(t_eme)
                for d = 1:length(sex_dep)
                    scedef.ID(sce) = sce;
                    scedef.time_app(sce) = time_app(a);
                    scedef.eff_size(sce) = eff_size(b);
                    scedef.t_eme(sce) = t_eme(c);
                    scedef.sex_dep(sce) = sex_dep(d);
            	    sce = sce + 1;
                end
            end
        end
    end

% Number of scenarios
nof_sce = size(scedef.ID,1);
sim_set.nof_sce = nof_sce;

% Save scenario definition
writetable(scedef,[dir.gensetdir '/scedef.csv'])

% Indicate scenarios to simulate
sce_mod = 1;
%sce_mod = [1,10,19,28,37,46,55];
    % 1 for control

loc_mod = 3;
%loc_mod = 1:3;

for l = loc_mod
    % Location parameters
    locparams = locations(l);

    % CREATE SUB-DIRECTORIES
    % Location Directory
    dir.locdir = [dir.workdir dir.dirSep char(locparams.locname)];
    system(['mkdir ' dir.locdir]);

    % Save location parameters
    writetable(locparams,[dir.locdir dir.dirSep 'locparams_' char(locparams.locname) '.csv'])

    for sce = sce_mod
        sce_ini = sce(1);
        sce_end = sce(end);

        for num_trial = 3:3
            exitcode = chirosEme_execute(dir,scedef,l,locparams,sce_ini,sce_end,num_trial);
        end
    end
end

return
read_results
plot_control