function [exitcode] = chirosEme_execute(dir,scedef,l,locparams,sce_ini,sce_end,num_trial)

	exitcode = 0;
    
    trial_dir = [dir.locdir  dir.dirSep 't' num2str(num_trial,'%03d')];
    system(['mkdir ' trial_dir]);

    % TIME SETTINGS
    % Time step (days)
    dt = 1/24; % 1 hour
    sim_set.dt = dt;
    
    % Start time (LST, Local Standard Time)
    dateStart_str = "01-May-2021 00:00:00";
    sim_set.dateStart_str = dateStart_str;

    dateStart = datenum(dateStart_str,"dd-mmm-yyyy HH:MM:SS");
    sim_set.dateStart = dateStart;

    t_start = dateStart;
    sim_set.t_start = t_start; % MATLAB serial date number
    
    % End time (LST, Local Standard Time)
    dateEnd_str = "01-Jan-2026 00:00:00";
    sim_set.dateEnd_str = dateEnd_str;

    dateEnd = datenum(dateEnd_str,"dd-mmm-yyyy HH:MM:SS");
    sim_set.dateEnd = dateEnd;

    t_end = dateEnd;
    sim_set.t_end = t_end; % MATLAB serial date number
    
    % Starting year
    year_start = year(dateStart_str);
    sim_set.year_start = year_start;
    
    year_start_Bti = 2024;
    sim_set.year_start_Bti = year_start_Bti;

    % End year
    year_end = year(dateEnd_str);
    sim_set.year_end = year_end;
    
    % Number of years
    nof_years = max(1,floor((t_end-t_start)/365)) + 1;
    sim_set.nof_years = nof_years;

    % Maximum number of time levels
    n_max = max(1,24*365*nof_years);
    sim_set.n_max = n_max;

    % DEFINITION OF Bti APPLICATION SEASON
    % First day of the Bti application season
    doy_Bti_start = 74; % 15th of March
    sim_set.doy_Bti_start = doy_Bti_start;

    % Last day of the Bti application season
    doy_Bti_end = 227; % 15th of August
    sim_set.doy_Bti_end = doy_Bti_end;

    % GENERAL SETTINGS
    % Number of life stages
    nof_stages = 6; % eggs + four larval instars + midges
    
    % Maximum number of age groups per life stage
    ag_max(1,1:2) = 1200;               % eggs
    ag_max(2:3,1:2) = 1200;             % 1st and 2nd instar larvae
    ag_max(4:5,1:2) = 4800;             % 3rd and 4th instar larvae
    ag_max(nof_stages,1) = 1200;        % male midges
    ag_max(nof_stages,2) = 1200;        % female midges
    
    % Age group areal density tolerance
    % (below this number, an age group is considered unpopulated)
    tol = 1e-4;
    
    % Ages of individuals in each age group
    age = cell(nof_stages,2);
    for s = 1:nof_stages
        age{s,1}(1:ag_max(s,1)) = (0:ag_max(s,1)-1)*dt;
        age{s,2}(1:ag_max(s,2)) = (0:ag_max(s,2)-1)*dt;
    end
    
    % SOURCES OF INDIVIDUALS IN EACH LIFE STAGE
    % Source 1: Life stage that supplies new individuals
    src1 = [6 1 2 3 4 5]'; % Both sexes
    % Source 2: Sex that supplies new male individuals
    src2(1:nof_stages,1) = [2 1 1 1 1 1]';
    % Source 2: Sex that supplies new female individuals
    src2(1:nof_stages,2) = [2 2 2 2 2 2]';
    
    %% MODEL PARAMETERS
    % PARAMETERS: GDD-BASED DEVELOPMENT MODEL
    % Base temperature for the GDD accumulation model (°C)
    T_base = 7.2;
    
    % Mean life stage requirements
    GDD_req(1,1) = 10; % male eggs
    GDD_req(1,2) = 10; % female eggs
    for s = 2:nof_stages-1
        GDD_req(s,1) = 70; % male larvae
        GDD_req(s,2) = 80; % female larvae
    end
    
    % Standard deviation factors
    fsigma(1) = 0.1; % males
    fsigma(2) = 0.1; % females   

    % MIDGE LIFESPAN MODEL - DISTRIBUTION PARAMETERS
    % Males
    age_mu(1) = 4;                   % days
    age_sd(1) = fsigma(1)*age_mu(1); % days
    % Females
    age_mu(2) = 6;                   % days
    age_sd(2) = fsigma(1)*age_mu(1); % days

    % EGG DEPOSITION MODEL PARAMETERS
    % Maximum unit egg laying rate (day^-1)
    switch(num_trial)
        case(1)
            r_egglay_max = 0.55;
        case(2)
            r_egglay_max = 0.55*1.2;
        case(3)
            r_egglay_max = 0.55*1.5;
    end

    % Half saturation constant (midge gross productivity curve)
    P_half = 0.15;
    
    % Fraction of eggs that produce male larvae
    P_eggs_m = 0.5;
    
    % Fraction of eggs that produce female larvae
    P_eggs_f = 1 - P_eggs_m;

    % Bti APPLICATION MODEL PARAMETERS
    P_threshold = 0.3;

    % Duration of the effects of a Bti application (days)
    Bti_duration = 1;
    
    % Minimum interval between consecutive Bti applications (days)
    Bti_interval = 14;

    % LOSS RATE PARAMETERS
    % Eggs
    l_egg_ref = 0.02;
    l_ref(1,1) = l_egg_ref;
    l_ref(1,2) = l_egg_ref;
    
    % Larvae
    switch(l)
        case(1)
            l_lrv_ref = 0.0014;
        case(2)
            l_lrv_ref = 0.00177;
        case(3)
            switch(num_trial)
                case(1)
                    l_lrv_ref = 0.00963;
                case(2)
                    l_lrv_ref = 0.00963*1.85;
                case(3)
                    l_lrv_ref = 0.00963*2.88;
            end
    end

    % Percent reduction of reference unit loss rate for subsequent life stages
    l_red_rate = 0.2;
    for s = 2:nof_stages-1
        l_ref(s,1) = (1-l_red_rate)^(s-2)*l_lrv_ref;
        l_ref(s,2) = (1-l_red_rate)^(s-2)*l_lrv_ref;
    end
    
    % Bti-related unit loss rate of newly hatched larvae
    l_lrv_Bti = 0.1;
    
    % EFFECTS OF Bti
    % Eggs (unaffected by Bti)
    l_Bti(1,1) = 0;
    l_Bti(1,2) = 0;

    % Larvae (affected by Bti)
    % Percent reduction of Bti-related unit loss rate for subsequent life stages
    l_Bti_red_rate = 0.2;
    for s = 2:nof_stages-1
        l_Bti(s,1) = (1-l_Bti_red_rate)^(s-2)*l_lrv_Bti;
        l_Bti(s,2) = (1-l_Bti_red_rate)^(s-2)*l_lrv_Bti;
    end

    % DIAPAUSE-RELATED LOSS RATES
    % Non-diapausal stages
    % Eggs
    l_diap(1,1) = 1;
    l_diap(1,2) = 1;
    % First instar 
    l_diap(2,1) = 1;
    l_diap(2,2) = 1;
    % Second instar
    l_diap(3,1) = 1;
    l_diap(3,2) = 1;
    
    % Diapausal stages
    % Percent reduction of unit loss rate during diapause of late-instar
    % larvae
    switch(l)
        case(1)
            l_red_diap = 0.5;
        case(2)
            l_red_diap = 0.5;
        case(3)
            l_red_diap = 0.5;
    end
    % Third instar
    l_diap(4,1) = l_red_diap*l_ref(4,1);
    l_diap(4,2) = l_red_diap*l_ref(4,2);
    % Fourth instar
    l_diap(5,1) = l_red_diap*l_ref(5,1);
    l_diap(5,2) = l_red_diap*l_ref(5,1);

    %% MODEL IMPLEMENTATION PARAMETERS
    % Initialise multiplicative factors
    k1 = nan(nof_stages,2);
    k2 = nan(nof_stages,2);
    
    % Multiplicative factors: Eggs
    k1(1,1) = P_eggs_m;
    k1(1,2) = P_eggs_f;
    
    k2(1,1) = 1;
    k2(1,2) = 1;
    
    % Multiplicative factors: Larvae
    for s = 2:nof_stages-1
        k1(s,1) = 1;
        k1(s,2) = 1;
        k2(s,1) = 1;
        k2(s,2) = 1;
    end

    % Multiplicative factors: Midges
    k1(nof_stages,1) = 1;
    k1(nof_stages,2) = 1;
    k2(nof_stages,1) = 0;
    k2(nof_stages,2) = 0;

    % RUN SIMULATION
    for sce = sce_ini:sce_end
        % Scenario Directory
        scedir = [trial_dir dir.dirSep 'sce_' num2str(sce,'%03d')];
		system(['mkdir ' scedir]);

        % Results Directory (csv)
	    csvdir = [scedir dir.dirSep 'csv'];
	    system(['mkdir ' csvdir]);

        % SCENARIO MODIFIERS
        % A) Bti APPLICATION TIMING
        if(scedef.time_app(sce) == "I") 
            % If Bti applications follow a fixed calendar, then...
            switch(l)
                case(1)
                    Bti_dates = ["12-Apr-2024 12:00:00";... % year 4
                                 "11-Jun-2024 12:00:00";...
                                 "16-Jul-2024 12:00:00";...
                                 "12-Apr-2025 12:00:00";... % year 5
                                 "11-Jun-2025 12:00:00";...
                                 "16-Jul-2025 12:00:00";...
                                 "12-Aug-2035 12:00:00"];
                case(2)
                    Bti_dates = ["28-Mar-2024 12:00:00";... % year 4
                                 "06-Jun-2024 12:00:00";...
                                 "11-Jul-2024 12:00:00";...
                                 "28-Mar-2025 12:00:00";... % year 5
                                 "06-Jun-2025 12:00:00";...
                                 "11-Jul-2025 12:00:00";...
                                 "12-Aug-2035 12:00:00"];
                case(3)
                    Bti_dates = ["14-Mar-2024 12:00:00";... % year 4
                                 "13-May-2024 12:00:00";...
                                 "16-Jun-2024 12:00:00";...
                                 "17-Jul-2024 12:00:00";...
                                 "14-Mar-2025 12:00:00";... % year 5
                                 "13-May-2025 12:00:00";...
                                 "16-Jun-2025 12:00:00";...
                                 "17-Jul-2025 12:00:00";...
                                 "12-Aug-2035 12:00:00"];
                case(4)
                    Bti_dates = ["28-Apr-2023 12:00:00";... % year 3
                                 "14-May-2023 12:00:00";...
                                 "10-Jun-2023 12:00:00";...
                                 "01-Jul-2023 12:00:00";...
                                 "28-Jul-2023 12:00:00";...
                                 "12-Aug-2023 12:00:00";...
                                 "28-Apr-2024 12:00:00";... % year 4
                                 "14-May-2024 12:00:00";...
                                 "10-Jun-2024 12:00:00";...
                                 "01-Jul-2024 12:00:00";...
                                 "28-Jul-2024 12:00:00";...
                                 "12-Aug-2024 12:00:00";...
                                 "28-Apr-2025 12:00:00";... % year 5
                                 "14-May-2025 12:00:00";...
                                 "10-Jun-2025 12:00:00";...
                                 "01-Jul-2025 12:00:00";...
                                 "28-Jul-2025 12:00:00";...
                                 "12-Aug-2025 12:00:00";...
                                 "12-Aug-2035 12:00:00"];
            end
            for d = 1:length(Bti_dates)
                t_add_Bti(d) = datenum(Bti_dates(d),"dd-mmm-yyyy HH:MM:SS");
            end
        else
            Bti_dates = [];
        end
        
        % B) EFFECT SIZE (CHANGE IN MORTALITY RATE OF NEWLY HATCHED LARVAE)
        switch(scedef.eff_size(sce))
            case("N/A")
                effSz = 1; % Not Applicable = no change in mortality rate
            case("m")
                % 10% of increase in mortality of newly hatched larvae
                effSz = 1.1;
            case("i")
                % 50% of increase in mortality of newly hatched larvae
                effSz = 1.5;
            case("s")
                % 90% of increase in mortality of newly hatched larvae
                effSz = 1.9;
        end
        
        % C) TIME OF EMERGENCE (CHANGE IN REQUIRED GDDs)
        switch(scedef.t_eme(sce))
            case("N/A")
                fGDD = 1; % Not Applicable = no change in GDD requirements
            case("A")
                 % Anticipated emergence
                %fGDD = 0.90;  % Required GDDs reduced by 10%
                fGDD = 0.80;  % Required GDDs reduced by 20%
            case("D")
                % Delayed emergence
                % fGDD = 1.1; % Required GDDs increased by 10%
                fGDD = 1.2; % Required GDDs increased by 20%
            case("NE")
                % No effect on timing of emergence
                fGDD = 1; % No change in GDD requirements
        end

        % D) SEX DEPENDENCE OF EFFECT STRENGTH
        % Modifies both effect size and GDDs required to emerge
        switch(scedef.sex_dep(sce))
            case("N/A")
                fM = 1;
                fF = 1;
            case("SM")
                fM = 1.3;
                fF = 1;
            case("SF")
                fM = 1;
                fF = 1.3;
            case("ND")
                fM = 1;
                fF = 1;
        end
    
        % INITIALISE VARIABLES
        initialise
    
        % INITIAL CONDITIONS
        initial_conditions
        
        tic
        % EVOLVE MODEL
        while t_instant(n) + dt < t_end

            % TEMPORAL COORDINATES
            % Maximum time level
            n_max = max(n,n_max);
        
            % Instant time (days)
            t_instant(n+1) = t_instant(n) + dt;
            
            % Day of year (day number)
	        doy_now(n+1) = day(datetime(datestr(t_instant(n+1))),'dayofyear');
            
	        % MATLAB Serial Date Number in Universal Time (UT1)
	        datenum_UT1(n+1) = t_instant(n+1) - locparams.t_zone/24;
            
            % ENVIRONMENTAL CUES
            % Duration of sunlight (hours)
            sunh(n+1) = sunlight(locparams.lat,locparams.long,datenum_UT1(n+1));
	        
	        % Air temperature (°C)
	        T_air(n+1) = airtemp(locparams.T_amp,locparams.omega,locparams.delta,locparams.T_0,t_instant(n+1));

            % Flag: Diapause conditions met (1 = yes; 0 = no)
            if((T_air(n) + T_air(n+1))/2 < locparams.T_diap || (sunh(n) + sunh(n+1))/2 < locparams.daylength_diap)
                diap(n+1) = 1;
            else
                diap(n+1) = 0;
            end
        
            % AREAL DENSITIES AND LIFESTAGE SUBTOTALS (individuals/m^2)
            %areal_density{life_stage,sex}(age_group,time_level)
            % All life stages
            for s = 1:nof_stages
                % First age group
	            ad{s,1}(1,n+1) = sum(k1(s,1)*p{src1(s),src2(s,1)}(:,n)*dt.*ad{src1(s),src2(s,1)}(:,n),'omitnan'); % males
	            ad{s,2}(1,n+1) = sum(k1(s,2)*p{src1(s),src2(s,2)}(:,n)*dt.*ad{src1(s),src2(s,2)}(:,n),'omitnan'); % females
            
                % All other age groups
                ad{s,1}(2:ag_max(s,1),n+1) =  (1 - (k2(s,1)*p{s,1}(1:ag_max(s,1)-1,n) + l{s,1}(1:ag_max(s,1)-1,n))*dt).*ad{s,1}(1:ag_max(s,1)-1,n); % males
                ad{s,2}(2:ag_max(s,2),n+1) =  (1 - (k2(s,2)*p{s,2}(1:ag_max(s,2)-1,n) + l{s,2}(1:ag_max(s,2)-1,n))*dt).*ad{s,2}(1:ag_max(s,2)-1,n); % females
                
                % Life stage subtotals (all age groups)
                ad_subtot{s,1}(n+1) = sum(ad{s,1}(:,n+1),'omitnan');             % males
                ad_subtot{s,2}(n+1) = sum(ad{s,2}(:,n+1),'omitnan');             % females
                ad_subtot{s,3}(n+1) = ad_subtot{s,1}(n+1) + ad_subtot{s,2}(n+1); % both
            end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            % AREAL DENSITIES OF TOTAL LARVAE
            ad_lrv_tot{1,1}(n+1) = 0; % males
            ad_lrv_tot{1,2}(n+1) = 0; % females
            for s = 2:nof_stages-1
                ad_lrv_tot{1,1}(n+1) = ad_lrv_tot{1,1}(n+1) + ad_subtot{s,1}(n+1); % males
                ad_lrv_tot{1,2}(n+1) = ad_lrv_tot{1,2}(n+1) + ad_subtot{s,2}(n+1); % females
            end
            ad_lrv_tot{1,3}(n+1) = ad_lrv_tot{1,1}(n+1) + ad_lrv_tot{1,2}(n+1); % both
        
            % LIFE STAGE RATIOS
            % LARVAE IN EACH INSTAR TO TOTAL LARVAE
            for s = 2:nof_stages-1
                P_lrv{s-1,1}(n+1) = ad_subtot{s,1}(n+1)/ad_lrv_tot{1,1}(n+1); % males
                P_lrv{s-1,2}(n+1) = ad_subtot{s,2}(n+1)/ad_lrv_tot{1,2}(n+1); % females
                P_lrv{s-1,3}(n+1) = ad_subtot{s,3}(n+1)/ad_lrv_tot{1,3}(n+1); % both
            end
        
            % LARVAE IN EACH LEVEL OF MATURITY TO TOTAL LARVAE
            % Proportion of young larvae
            P_lrv_yng{1,1}(n+1) = (ad_subtot{2,1}(n+1) + ad_subtot{3,1}(n+1))/ad_lrv_tot{1,1}(n+1);	% males
            P_lrv_yng{1,2}(n+1) = (ad_subtot{2,2}(n+1) + ad_subtot{3,2}(n+1))/ad_lrv_tot{1,2}(n+1);	% females
            P_lrv_yng{1,3}(n+1) = (ad_subtot{2,3}(n+1) + ad_subtot{3,3}(n+1))/ad_lrv_tot{1,3}(n+1);	% both
            % Proportion of mature larvae
            P_lrv_mat{1,1}(n+1) = (ad_subtot{4,1}(n+1) + ad_subtot{5,1}(n+1))/ad_lrv_tot{1,1}(n+1);	% males
            P_lrv_mat{1,2}(n+1) = (ad_subtot{4,2}(n+1) + ad_subtot{5,2}(n+1))/ad_lrv_tot{1,2}(n+1);	% females
            P_lrv_mat{1,3}(n+1) = (ad_subtot{4,3}(n+1) + ad_subtot{5,3}(n+1))/ad_lrv_tot{1,3}(n+1);	% both
        
            % SEX RATIOS: ALL LIFE STAGES
            for s = 1:nof_stages
                % Male individuals in life stage 's' to total individuals in the same life stage
                P_m2a{s,1}(n+1) = ad_subtot{s,1}(n+1)/ad_subtot{s,3}(n+1);
                % Female individuals in life stage 's' to total individuals in the same life stage
                P_f2a{s,1}(n+1) = ad_subtot{s,2}(n+1)/ad_subtot{s,3}(n+1);
                % Male individuals in life stage 's' to female individuals in the same life stage
                P_m2f{s,1}(n+1) = ad_subtot{s,1}(n+1)/ad_subtot{s,2}(n+1);
            end

            % SEX RATIOS: LARVAE ONLY
            % Young male larvae to total young larvae
            P_m2a_lrv_yng(n+1) = (ad_subtot{2,1}(n+1) + ad_subtot{3,1}(n+1))/(ad_subtot{2,3}(n+1) + ad_subtot{3,3}(n+1));
            % Mature male larvae to total mature larvae
            P_m2a_lrv_mat(n+1) = (ad_subtot{4,1}(n+1) + ad_subtot{5,1}(n+1))/(ad_subtot{4,3}(n+1) + ad_subtot{5,3}(n+1));
            % All male larvae to total larvae
            P_m2a_lrv_all(n+1) = ad_lrv_tot{1,1}(n+1)/ad_lrv_tot{1,3}(n+1);
            
            % Young female larvae to total young larvae
            P_f2a_lrv_yng(n+1) = (ad_subtot{2,2}(n+1) + ad_subtot{3,2}(n+1))/(ad_subtot{2,3}(n+1) + ad_subtot{3,3}(n+1));
            % Mature female larvae to total mature larvae
            P_f2a_lrv_mat(n+1) = (ad_subtot{4,2}(n+1) + ad_subtot{5,2}(n+1))/(ad_subtot{4,3}(n+1) + ad_subtot{5,3}(n+1));
            % All female larvae to total larvae
            P_f2a_lrv_all(n+1) = ad_lrv_tot{1,2}(n+1)/ad_lrv_tot{1,3}(n+1);
            
            % Young male larvae to young female larvae
            P_m2f_lrv_yng(n+1) = (ad_subtot{2,1}(n+1) + ad_subtot{3,1}(n+1))/(ad_subtot{2,2}(n+1) + ad_subtot{3,2}(n+1));
            % Mature male larvae to mature female larvae
            P_m2f_lrv_mat(n+1) = (ad_subtot{4,1}(n+1) + ad_subtot{5,1}(n+1))/(ad_subtot{4,2}(n+1) + ad_subtot{5,2}(n+1));
            % All male larvae to all female larvae
            P_m2f_lrv_all(n+1) = ad_lrv_tot{1,1}(n+1)/ad_lrv_tot{1,2}(n+1);
            
            %---------------------%
            % AUXILIARY VARIABLES %
            %---------------------%
            % Bti APPLICATION VARIABLES
            % Flag: Bti application season (1 = on; 0 = off)
            if(doy_now(n+1) >= doy_Bti_start && doy_now(n+1) <= doy_Bti_end && year(t_instant(n+1)) >= year_start_Bti)
                Bti_season(n+1) = 1;
            else
                Bti_season(n+1) = 0;
            end

            % Flag: Bti application effects (1 = yes; 0 = no)
            % By default, equal to the value of the previous time step
            Bti_apply(n+1) = Bti_apply(n);
            
            % Counter: Index number of the next Bti application
            % By default, equal to the value of the previous time step
            Bti_app_no(n+1) = Bti_app_no(n);

            % Is this the "Control" scenario?
            if(scedef.time_app(sce) ~= "N/A")
                % If this is NOT the "Control" scenario, then...
                
                % Are the effects of a Bti application active?
                if(Bti_apply(n+1))
                    % YES, the effects of a Bti application are active...

                    % But should these effects remain active?
                    if(t_instant(n+1) > t_app_end)
                        % If the effects of the last Bti application have been active for 24 hours or more, then...
                        % dectivate application of B.t.i.
                        Bti_apply(n+1) = 0;
                    end

                else
                    % NO, the effects of a Bti application are NOT active...

                    % DECISION: Apply Bti?
                    if(scedef.time_app(sce) == "I")
                        % If the application of Bti follows a fixed calendar, then ...
                        
                        if(t_instant(n+1) >= t_add_Bti(Bti_app_no(n+1)))
                            % If the new time is the first one to match or
                            % fall right after the scheduled time of the next
                            % application, then...

                            % Activate effects of the application of Bti
                            Bti_apply(n+1) = 1;

                            % End time of the effects of this application
                            t_app_end = t_instant(n+1) + Bti_duration;

                            % Increase the index number of the next application
                            Bti_app_no(n+1) = Bti_app_no(n+1) + 1;

                            if(Bti_app_no(n+1) > length(t_add_Bti))
                                % If all scheduled applications have happened, then ... 
                                % no further applications will happen
                                Bti_app_no(n+1) = NaN;
                            end

                        end

                    elseif(Bti_season(n+1) && ...
                                ((scedef.time_app(sce) == "Y" && ...
                                  P_lrv_yng{1,3}(n-1) < P_lrv_yng{1,3}(n) && ...
                                  P_lrv_yng{1,3}(n+1) < P_lrv_yng{1,3}(n) && ...
                                  P_lrv_yng{1,3}(n+1) > P_lrv_mat{1,3}(n+1) && ...
                                  P_lrv_yng{1,3}(n+1) > P_threshold) || ...
                                 (scedef.time_app(sce) == "X" && ...
                                  P_lrv_yng{1,3}(n-1) + P_lrv_mat{1,3}(n-1) < P_lrv_yng{1,3}(n) + P_lrv_mat{1,3}(n) && ...
                                  P_lrv_yng{1,3}(n+1) + P_lrv_mat{1,3}(n+1) < P_lrv_yng{1,3}(n) + P_lrv_mat{1,3}(n)) || ...
                                 (scedef.time_app(sce) == "M" && ...
                                  P_lrv_mat{1,3}(n-1) < P_lrv_mat{1,3}(n) && ...
                                  P_lrv_mat{1,3}(n+1) < P_lrv_mat{1,3}(n) && ...
                                  P_lrv_mat{1,3}(n+1) > P_lrv_yng{1,3}(n+1) && ...
                                  P_lrv_mat{1,3}(n+1) > P_threshold)))
                                
                                % If the application of Bti happens under certain scenario-specific conditions and ...
                                % these scenario-specific conditions are met, then...

                                if(Bti_app_no(n+1) == 1)% Activate effects of the application of Bti
                                    Bti_apply(n+1) = 1;
                                    
                                    % Update list of Bti application dates
                                    Bti_dates = [Bti_dates; string(datestr(t_instant(n+1)))];
                                    
                                    % End time of the effects of this application
                                    t_app_end = t_instant(n+1) + Bti_duration;
            
                                    % Increase the index number of the next application
                                    Bti_app_no(n+1) = Bti_app_no(n+1) + 1;
        
                                elseif(t_instant(n+1) - datenum(Bti_dates(end)) > Bti_interval)
                                    Bti_apply(n+1) = 1;
        
                                    % Update list of Bti application dates
                                    Bti_dates = [Bti_dates; string(datestr(t_instant(n+1)))];
                                    
                                    % End time of the effects of this application
                                    t_app_end = t_instant(n+1) + Bti_duration;
            
                                    % Increase the index number of the next application
                                    Bti_app_no(n+1) = Bti_app_no(n+1) + 1;
        
                                end
                    end
                end
            end

            % GDD-BASED DEVELOPMENT AND EVOLUTION MODEL
            % GDDs accumulated over the last time step
	        dGDD(n+1) = max(0,((T_air(n) + T_air(n+1))/2 - T_base)*dt);

            % GDDs accumulated by each age group since initiating the life stage
            for s = 1:nof_stages-1  
                % First age group
                GDD{s,1}(1,n+1) = dGDD(n+1); % males
                GDD{s,2}(1,n+1) = dGDD(n+1); % females
            
                % All other age groups
                GDD{s,1}(2:ag_max(s,1),n+1) = GDD{s,1}(1:ag_max(s,1)-1,n) + dGDD(n+1);
                GDD{s,2}(2:ag_max(s,2),n+1) = GDD{s,2}(1:ag_max(s,2)-1,n) + dGDD(n+1);
            end

            % Distribution parameters of GDD-based model
            for s = 1:nof_stages-1    
                % Mean
	            GDD_mu{s,1}(n+1) = (Bti_apply(n+1)*fGDD*fM + (1 - Bti_apply(n+1)))*GDD_req(s,1); % males
                GDD_mu{s,2}(n+1) = (Bti_apply(n+1)*fGDD*fF + (1 - Bti_apply(n+1)))*GDD_req(s,2); % females
            
                % Standard deviation
                GDD_sd{s,1}(n+1) = fsigma(1)*GDD_mu{s,1}(n+1); % males
                GDD_sd{s,2}(n+1) = fsigma(2)*GDD_mu{s,2}(n+1); % females
            end
            
            % GROSS PRODUCTIVITIES
            % Eggs and larvae
            for s = 1:nof_stages-1
                p{s,1}(:,n+1) = 0 + (1 - diap(n+1))*normcdf(GDD{s,1}(:,n+1),GDD_mu{s,1}(n+1),GDD_sd{s,1}(n+1))';
                p{s,2}(:,n+1) = 0 + (1 - diap(n+1))*normcdf(GDD{s,2}(:,n+1),GDD_mu{s,2}(n+1),GDD_sd{s,2}(n+1))';
            end
            % Midges
            p{nof_stages,1}(:,n+1) = 0;
            p{nof_stages,2}(:,n+1) = r_egglay_max*P_m2f{nof_stages,1}(n+1)/(P_half + P_m2f{nof_stages,1}(n+1));
            
            % LOSS RATES
            % Eggs and larvae
            for s = 1:nof_stages-1
                % First age group
                l{s,1}(1,n+1) = (1 - diap(n+1))*((1 - Bti_apply(n+1))*l_ref(s,1) + Bti_apply(n+1)*effSz*fM*l_Bti(s,1)) + ...
                                (    diap(n+1))*l_diap(s,1); % males
                l{s,2}(1,n+1) = (1 - diap(n+1))*((1 - Bti_apply(n+1))*l_ref(s,2) + Bti_apply(n+1)*effSz*fF*l_Bti(s,2)) + ...
                                (    diap(n+1))*l_diap(s,2); % females
                % All other age groups
                l{s,1}(2:ag_max(s,1),n+1) = [0.999.^((2:ag_max(s,1))-2)*l{s,1}(1,n+1)]'; % males
                l{s,2}(2:ag_max(s,2),n+1) = [0.999.^((2:ag_max(s,2))-2)*l{s,2}(1,n+1)]'; % females
            end
            % Midges
            l{nof_stages,1}(1:ag_max(nof_stages,1),n+1) = normcdf(age{nof_stages,1},age_mu(1),age_sd(1));
            l{nof_stages,2}(1:ag_max(nof_stages,2),n+1) = normcdf(age{nof_stages,2},age_mu(2),age_sd(2));
        
            % Make NaN all rates for unpopulated age groups within each stage
            for s = 1:nof_stages
                p{s,1}(isnan(ad{s,1}(:,n+1)) | ad{s,1}(:,n+1) < 0,n+1) = NaN;
                p{s,2}(isnan(ad{s,2}(:,n+1)) | ad{s,2}(:,n+1) < 0,n+1) = NaN;
                l{s,1}(isnan(ad{s,1}(:,n+1)) | ad{s,1}(:,n+1) < 0,n+1) = NaN;
                l{s,2}(isnan(ad{s,2}(:,n+1)) | ad{s,2}(:,n+1) < 0,n+1) = NaN;
            end
        
%             for s = 1:nof_stages
%                 if(ad{s,1}(ag_max(s,1),n+1) > tol)
%                     ag_max(s,1) = ag_max(s,1) + 1;
%                     age{s,1}(1:ag_max(s,1)) = (0:ag_max(s,1)-1)*dt;
%                 end
%                 if(ad{s,2}(ag_max(s,2),n+1) > tol)
%                     ag_max(s,2) = ag_max(s,2) + 1;
%                     age{s,2}(1:ag_max(s,2)) = (0:ag_max(s,2)-1)*dt;
%                 end
%             end
        
            % Update time level
            n = n + 1;
        end
        toc
        
        % SAVE SIMULATION RESULTS
        dirSep = dir.dirSep;
        savesim
        
        exitcode = 1;

        % CLEAR VARIABLES
        clear t_instant doy_now datenum_UT1 sunh T_air diap ...
              ad ad_subtot ad_lrv_tot p l...
              P_lrv P_lrv_yng P_lrv_mat ...
              P_m2a P_f2a P_m2f ...
              P_m2a_lrv_yng P_m2a_lrv_mat P_m2a_lrv_all ...
              dGDD GDD GDD_mu GDD_sd ...
              Bti_season Bti_apply Bti_app_no Bti_dates
    end
end