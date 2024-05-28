%% INITIALISE VARIABLES
% Initialise vectors and arrays
t_instant = nan(1,n_max);
doy_now = nan(1,n_max);
datenum_UT1 = nan(1,n_max);
sunh = nan(1,n_max);
T_air = nan(1,n_max);
dGDD = nan(1,n_max);
diap = nan(1,n_max);
Bti_season = nan(1,n_max);
Bti_apply = nan(1,n_max);
Bti_app_no = nan(1,n_max);

% GDD-based development model
GDD = cell(nof_stages-1,2);
GDD_mu = cell(nof_stages-1,2);
GDD_sd = cell(nof_stages-1,2);
for s = 1:nof_stages-1
    GDD{s,1}(1:ag_max(s,1),1:n_max) = NaN;
    GDD{s,2}(1:ag_max(s,2),1:n_max) = NaN;
    GDD_mu{s,1}(1:n_max) = NaN;
    GDD_mu{s,2}(1:n_max) = NaN;
    GDD_sd{s,1}(1:n_max) = NaN;
    GDD_sd{s,2}(1:n_max) = NaN;
end

% AREAL DENSITIES
ad = cell(nof_stages,2);
ad_subtot = cell(nof_stages,3);
for s = 1:nof_stages
    ad{s,1}(1:ag_max(s,1),1:n_max) = NaN;
    ad{s,2}(1:ag_max(s,2),1:n_max) = NaN;
    ad_subtot{s,1}(1:n_max) = NaN;
    ad_subtot{s,2}(1:n_max) = NaN;
    ad_subtot{s,3}(1:n_max) = NaN;
end
ad_lrv_tot = cell(nof_stages-2,3);
for s = 2:nof_stages-2
    ad_lrv_tot{s,1}(1:n_max) = NaN;
    ad_lrv_tot{s,2}(1:n_max) = NaN;
    ad_lrv_tot{s,3}(1:n_max) = NaN;
end

% POPULATION RATIOS
P_lrv = cell(nof_stages-2,3);
for s = 2:nof_stages-1
    P_lrv{s-1,1}(1:n_max) = NaN;
    P_lrv{s-1,2}(1:n_max) = NaN;
    P_lrv{s-1,3}(1:n_max) = NaN;
end

P_lrv_yng = cell(1,3);
P_lrv_yng{1,1}(1:n_max) = NaN;
P_lrv_yng{1,2}(1:n_max) = NaN;
P_lrv_yng{1,3}(1:n_max) = NaN;

P_lrv_mat = cell(1,3);
P_lrv_mat{1,1}(1:n_max) = NaN;
P_lrv_mat{1,2}(1:n_max) = NaN;
P_lrv_mat{1,3}(1:n_max) = NaN;

P_m2a = cell(nof_stages,1);
P_f2a = cell(nof_stages,1);
P_m2f = cell(nof_stages,1);
for s = 1:nof_stages
    P_m2a{s,1}(1:n_max) = NaN;
    P_f2a{s,1}(1:n_max) = NaN;
    P_m2f{s,1}(1:n_max) = NaN;
end

P_m2a_lrv_yng = nan(1,n_max);
P_m2a_lrv_mat = nan(1,n_max);
P_m2a_lrv_all = nan(1,n_max);

P_f2a_lrv_yng = nan(1,n_max);
P_f2a_lrv_mat = nan(1,n_max);
P_f2a_lrv_all = nan(1,n_max);

P_m2f_lrv_yng = nan(1,n_max);
P_m2f_lrv_mat = nan(1,n_max);
P_m2f_lrv_all = nan(1,n_max);

% Rates
p = cell(nof_stages,2);
l = cell(nof_stages,2);
for s = 1:nof_stages
    p{s,1}(1:ag_max(s,1),1:n_max) = NaN;
    p{s,2}(1:ag_max(s,2),1:n_max) = NaN;
    l{s,1}(1:ag_max(s,1),1:n_max) = NaN;
    l{s,2}(1:ag_max(s,2),1:n_max) = NaN;
end