%% INITIAL CONDITIONS

% TEMPORAL COORDINATES
% Initial time level
n = 1;

% Instant time (days)
t_instant(n) = t_start;

% Day of year (day number)
doy_now(n) = day(datetime(datestr(t_instant(n))),'dayofyear');

% MATLAB Serial Date Number in Universal Time (UT1)
datenum_UT1(n) = t_instant(n) - locparams.t_zone/24;

% ENVIRONMENTAL CUES
% Duration of sunlight (hours)
sunh(n) = sunlight(locparams.lat,locparams.long,datenum_UT1(n));
	
% Air temperature (°C)
T_air(n) = airtemp(locparams.T_amp,locparams.omega,locparams.delta,locparams.T_0,t_instant(n));

% Flag: Diapause conditions met (1 = yes; 0 = no)
if(T_air(n) < locparams.T_diap || sunh(n) < locparams.daylength_diap)
    diap(n) = 1;
else
    diap(n) = 0;
end

% AREAL DENSITIES (individuals/m^2)
% Males
% Areal density of newly deposed male chironomid eggs (eggs/m^2)
ad{1,1}(1,n) = 1.25e5;
% Areal densities of newly eclosed male chironomid larvae (larvae/m^2)
for s = 2:nof_stages-1
    ad{s,1}(1,n) = 0;
end
% Areal density of newly emerged male chironomid midges (midges/m^2)
ad{nof_stages,1}(1,n) = 0;

% Females
% Areal density of newly deposed female chironomid eggs (eggs/m^2)
ad{1,2}(1,n) = 1.25e5;
% Areal densities of newly eclosed female chironomid larvae (larvae/m^2)
for s = 2:nof_stages-1
    ad{s,2}(1,n) = 0;
end
% Areal density of newly emerged female chironomid midges (midges/m^2)
ad{nof_stages,2}(1,n) = 0;

% SUBTOTALS PER STAGE (all age groups within each life stage)
for s = 1:nof_stages
    ad_subtot{s,1}(n) = sum(ad{s,1}(:,n),'omitnan');            % males
    ad_subtot{s,2}(n) = sum(ad{s,2}(:,n),'omitnan');            % females
    ad_subtot{s,3}(n) = ad_subtot{s,1}(n) + ad_subtot{s,2}(n);  % both
end

% TOTAL NUMBER OF LARVAE (all age groups in all larval stages)
ad_lrv_tot{1,1}(n) = 0; % males
ad_lrv_tot{1,2}(n) = 0; % females
for s = 2:nof_stages-1
    ad_lrv_tot{1,1}(n) = ad_lrv_tot{1,1}(n) + ad_subtot{s,1}(n); % males
    ad_lrv_tot{1,2}(n) = ad_lrv_tot{1,2}(n) + ad_subtot{s,2}(n); % females
end
ad_lrv_tot{1,3}(n) = ad_lrv_tot{1,1}(n) + ad_lrv_tot{1,2}(n); % both

% LIFE STAGE RATIOS
% LARVAE IN EACH INSTAR TO TOTAL LARVAE
for s = 2:nof_stages-1
    P_lrv{s-1,1}(n) = ad_subtot{s,1}(n)/ad_lrv_tot{1,1}(n); % males
    P_lrv{s-1,2}(n) = ad_subtot{s,2}(n)/ad_lrv_tot{1,2}(n); % females
    P_lrv{s-1,3}(n) = ad_subtot{s,3}(n)/ad_lrv_tot{1,3}(n); % both
end

% LARVAE IN EACH LEVEL OF MATURITY TO TOTAL LARVAE
% Proportion of young larvae
P_lrv_yng{1,1}(n) = (ad_subtot{2,1}(n) + ad_subtot{3,1}(n))/ad_lrv_tot{1,1}(n);	% males
P_lrv_yng{1,2}(n) = (ad_subtot{2,2}(n) + ad_subtot{3,2}(n))/ad_lrv_tot{1,2}(n);	% females
P_lrv_yng{1,3}(n) = (ad_subtot{2,3}(n) + ad_subtot{3,3}(n))/ad_lrv_tot{1,3}(n);	% both
% Proportion of mature larvae
P_lrv_mat{1,1}(n) = (ad_subtot{4,1}(n) + ad_subtot{5,1}(n))/ad_lrv_tot{1,1}(n);	% males
P_lrv_mat{1,2}(n) = (ad_subtot{4,2}(n) + ad_subtot{5,2}(n))/ad_lrv_tot{1,2}(n);	% females
P_lrv_mat{1,3}(n) = (ad_subtot{4,3}(n) + ad_subtot{5,3}(n))/ad_lrv_tot{1,3}(n);	% both

% SEX RATIOS: ALL LIFE STAGES
for s = 1:nof_stages
    % Male individuals in life stage 's' to total individuals in the same life stage
    P_m2a{s,1}(n) = ad_subtot{s,1}(n)/ad_subtot{s,3}(n);
    % Female individuals in life stage 's' to total individuals in the same life stage
    P_f2a{s,1}(n) = ad_subtot{s,2}(n)/ad_subtot{s,3}(n);
    % Male individuals in life stage 's' to female individuals in the same life stage
    P_m2f{s,1}(n) = ad_subtot{s,1}(n)/ad_subtot{s,2}(n);
end

% SEX RATES: LARVAE ONLY
% Young male larvae to total young larvae
P_m2a_lrv_yng(n) = (ad_subtot{2,1}(n) + ad_subtot{3,1}(n))/(ad_subtot{2,3}(n) + ad_subtot{3,3}(n));
% Mature male larvae to total mature larvae
P_m2a_lrv_mat(n) = (ad_subtot{4,1}(n) + ad_subtot{5,1}(n))/(ad_subtot{4,3}(n) + ad_subtot{5,3}(n));
% All male larvae to total larvae
P_m2a_lrv_all(n) = ad_lrv_tot{1,1}(n)/ad_lrv_tot{1,3}(n);

% Young female larvae to total young larvae
P_f2a_lrv_yng(n) = (ad_subtot{2,2}(n) + ad_subtot{3,2}(n))/(ad_subtot{2,3}(n) + ad_subtot{3,3}(n));
% Mature female larvae to total mature larvae
P_f2a_lrv_mat(n) = (ad_subtot{4,2}(n) + ad_subtot{5,2}(n))/(ad_subtot{4,3}(n) + ad_subtot{5,3}(n));
% All female larvae to total larvae
P_f2a_lrv_all(n) = ad_lrv_tot{1,2}(n)/ad_lrv_tot{1,3}(n);

% Young male larvae to young female larvae
P_m2f_lrv_yng(n) = (ad_subtot{2,1}(n) + ad_subtot{3,1}(n))/(ad_subtot{2,2}(n) + ad_subtot{3,2}(n));
% Mature male larvae to mature female larvae
P_m2f_lrv_mat(n) = (ad_subtot{4,1}(n) + ad_subtot{5,1}(n))/(ad_subtot{4,2}(n) + ad_subtot{5,2}(n));
% All male larvae to all female larvae
P_m2f_lrv_all(n) = ad_lrv_tot{1,1}(n)/ad_lrv_tot{1,2}(n);

%---------------------%
% AUXILIARY VARIABLES %
%---------------------%
% Bti APPLICATION VARIABLES
% Flag: B.t.i. application season (1 = on; 0 = off)
Bti_season(n) = 0;

% Flag: Bti application effects (1 = yes; 0 = no)
Bti_apply(n) = 0;

% Counter: Index number of the next Bti application
Bti_app_no(n) = 1;

% GDD-BASED DEVELOPMENT AND EVOLUTION MODEL
% GDDs accumulated over the last time step
dGDD(n) = 0;

for s = 1:nof_stages-1
    % Accumulated GDDs since initiating the life stage
    GDD{s,1}(~isnan(ad{s,1}(:,n)) & ad{s,1}(:,n) > 0,n) = 0;
    GDD{s,2}(~isnan(ad{s,2}(:,n)) & ad{s,2}(:,n) > 0,n) = 0;

    % Distribution parameters of GDD-based model
    % Mean
	GDD_mu{s,1}(n) = (Bti_apply(n)*fGDD*effSz*fM + (1 - Bti_apply(n)))*GDD_req(s,1); % males
    GDD_mu{s,2}(n) = (Bti_apply(n)*fGDD*effSz*fF + (1 - Bti_apply(n)))*GDD_req(s,2); % females

    % Standard deviation
    GDD_sd{s,1}(n) = fsigma(1)*GDD_mu{s,1}(n); % males
    GDD_sd{s,2}(n) = fsigma(2)*GDD_mu{s,2}(n); % females
end

% GROSS PRODUCTIVITIES
% Eggs and larvae
for s = 1:nof_stages-1
    p{s,1}(:,n) = 0 + (1 - diap(n))*normcdf(GDD{s,1}(:,n),GDD_mu{s,1}(n),GDD_sd{s,1}(n));
    p{s,2}(:,n) = 0 + (1 - diap(n))*normcdf(GDD{s,2}(:,n),GDD_mu{s,2}(n),GDD_sd{s,2}(n));
end
% Midges
p{nof_stages,1}(~isnan(ad{nof_stages,1}(:,n)) & ad{nof_stages,1}(:,n) >= 0,n) = 0;
p{nof_stages,2}(~isnan(ad{nof_stages,1}(:,n)) & ad{nof_stages,1}(:,n) >= 0,n) = r_egglay_max*P_m2f{nof_stages,1}(n)/(P_half + P_m2f{nof_stages,1}(n));

% LOSS RATES
% Eggs and Larvae
for s = 1:nof_stages-1
    l{s,1}(1,n) = (1 - diap(n))*(l_ref(s,1) + Bti_apply(n)*effSz*fM*l_Bti(s,1)) + ...
                       diap(n)*(l_diap(s,1) + Bti_apply(n)*effSz*fM*l_Bti(s,1));
    l{s,2}(1,n) = (1 - diap(n))*(l_ref(s,2) + Bti_apply(n)*effSz*fF*l_Bti(s,2)) + ...
                       diap(n)*(l_diap(s,2) + Bti_apply(n)*effSz*fF*l_Bti(s,2));
    for ag = 2:ag_max(s,1)
        l{s,1}(ag,n) = 0.999^(ag-2)*l{s,1}(1,n);
    end
    for ag = 2:ag_max(s,2)
        l{s,2}(ag,n) = 0.999^(ag-2)*l{s,2}(1,n);
    end
end
% Midges
l{nof_stages,1}(1:ag_max(nof_stages,1),n) = normcdf(age{nof_stages,1},age_mu(1),age_sd(1));
l{nof_stages,2}(1:ag_max(nof_stages,2),n) = normcdf(age{nof_stages,2},age_mu(2),age_sd(2));

for s = 1:nof_stages
    p{s,1}(isnan(ad{s,1}(:,n)) | ad{s,1}(:,n) < 0,n) = NaN;
    p{s,2}(isnan(ad{s,2}(:,n)) | ad{s,2}(:,n) < 0,n) = NaN;
    l{s,1}(isnan(ad{s,1}(:,n)) | ad{s,1}(:,n) < 0,n) = NaN;
    l{s,2}(isnan(ad{s,2}(:,n)) | ad{s,2}(:,n) < 0,n) = NaN;
end

% End of application effects
t_app_end = NaN;