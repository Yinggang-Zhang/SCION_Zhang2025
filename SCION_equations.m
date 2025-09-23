function dy = SCION_equations(t,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%              110111010                                                                        %
%           111010-1-----101                                                                    %
%        1011111---------101111                                                                 %
%      11011------------------101         SCION: Spatial Continuous Integration                 %
%     111-----------------10011011        Earth Evolution Model                                 %
%    1--10---------------1111011111                                                             %
%    1---1011011---------1010110111       Coded by Benjamin J. W. Mills                         %
%    1---1011000111----------010011       email: b.mills@leeds.ac.uk                            %
%    1----1111011101----------10101                                                             %
%     1----1001111------------0111        Model equations file                                  %
%      1----1101-------------1101         contains flux and reservoir equations                 %
%        1--111----------------1                                                                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% setup dy array
dy = zeros(23,1);  

%%%%%%% set up global parameters
global stepnumber
global pars
global forcings
global workingstate
global gridstate
global INTERPSTACK
global sensanal
global sensparams
global Insolation_mode
global time_index
global insolation_lat_variations

%%%%%%%%%%%%% get variables from Y to make working easier
P_i = y(1) ;
P_s = y(22) ;

O = y(2) ;
A = y(3) ;
S = y(4) ;
G = y(5) ;
C = y(6) ;
PYR = y(7) ;
GYP = y(8) ;
% TEMP = y(9);
% CAL = y(10) ;
N_i = y(11) ;
N_s = y(23) ;

OSr = y(18) ;
SSr = y(20) ;
dSSr = y(21)/y(20) ;

%%%% geological time in Ma (S=-2 sets a present day run)
if pars.runcontrol == -2
    t_geol = 0 ;
else
    t_geol = t*(1e-6) ;
end

%%%%%%% calculate isotopic fractionation of reservoirs
delta_G = y(12)/y(5);
delta_C = y(13)/y(6);
delta_GYP  = y(15)/y(8);
delta_PYR  = y(14)/y(7);

%%%%%%% atmospheric fraction of total CO2, atfrac(A)
atfrac0 = 0.01614 ;
%%%%%%% constant
% atfrac = 0.01614 ;
%%%%%%% variable
atfrac = atfrac0 * (A/pars.A0) ;

%%%%%%%% calculations for pCO2, pO2
RCO2 = (A/pars.A0)*(atfrac/atfrac0) ;
CO2atm = RCO2*(280e-6) ;
CO2ppm = RCO2*280 ;

%%%%% mixing ratio of oxygen (not proportional to O reservoir)
mrO2 = ( O/pars.O0 )  /   ( (O/pars.O0)  + pars.copsek16 ) ;
%%%%% relative moles of oxygen 
RO2 =  O/pars.O0 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Interpolate forcings for this timestep   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% COPSE Reloaded forcing set
E_reloaded = interp1( forcings.t, forcings.E , t_geol ) ;
W_reloaded = interp1( forcings.t, forcings.W , t_geol ) ;
%%%% Additional forcings
GR_BA = interp1( forcings.GR_BA(:,1)./1e6 , forcings.GR_BA(:,2) , t_geol ) ;
newGA = interp1( forcings.newGA(:,1)./1e6 , forcings.newGA(:,2) , t_geol ) ;
D_combined_mid = interp1( forcings.D_force_x , forcings.D_force_mid, t_geol) ;
D_combined_min = interp1( forcings.D_force_x , forcings.D_force_min, t_geol) ;
D_combined_max = interp1( forcings.D_force_x , forcings.D_force_max, t_geol) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Choose forcing functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEGASS = 1 ;
% DEGASS = D_sbz_rift ;
% DEGASS = D_merdith ;
DEGASS = D_combined_mid ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = 1 ;
W = W_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVO = 1 ;
EVO = E_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CPLAND = 1 ;
% CPLAND = CP_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bforcing = interp1qr([-1000 -150 -100 0]',[0.75 0.75 1 1]',t_geol) ;
% Bforcing = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BAS_AREA = GR_BA ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAN_AREA = newGA ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PREPLANT = 1/4 ;
capdelS = 27   ;
capdelC_land = 27   ;
capdelC_marine = 35  ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% enhancement of inso-temperature factor
factor1 = 3;
%%%%%%%%%%


%%%% SHORELINE
SHORELINE = interp1qr(forcings.shoreline_time',forcings.shoreline_relative',t_geol) ;

%%%% bioturbation forcing
% f_biot = interp1qr([-1000 -525 -520 0]',[0 0 1 1]',t_geol);
% CB = interp1qr([0 1]',[1.2 1]',f_biot) ;
f_biot = interp1qr([-1000 -525 -520 0]',[1 1 1 1]',t_geol);
CB = interp1qr([0 1]',[1 1]',f_biot) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Sensitivity analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% all sensparams vary between [-1 +1]
if sensanal == 1
    
    %%%% Very degassing between upper and lower bounds
    if sensparams.randminusplus1 > 0
        DEGASS = (1 - sensparams.randminusplus1)*DEGASS + sensparams.randminusplus1*D_combined_max ;
    else
        DEGASS = (1 + sensparams.randminusplus1)*DEGASS - sensparams.randminusplus1*D_combined_min ;
    end
        
    %%%% simple +/- 20% variation
    BAS_AREA = BAS_AREA * (1 + 0.2*sensparams.randminusplus2) ;
    GRAN_AREA = GRAN_AREA * (1 + 0.2*sensparams.randminusplus3) ;
    
    %%%% preplant varies from 1/1 to 1/7
    PREPLANT = 1/ ( 4 + 3*sensparams.randminusplus4) ;
    
    %%%%
    capdelS = 30 + 5*sensparams.randminusplus5 ;
    capdelC_land = 25 + 5*sensparams.randminusplus6 ;
    capdelC_marine = 30 + 5*sensparams.randminusplus7 ;
    
    %%% 20% uncertain in tem-inso factor
    factor1 = factor1 * (1 + 0.2 * sensparams.randminusplus8) ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial fields from stack   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Fetch keyframe grids   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% find past and future keyframes
key_past_time = max ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) <= 0 ) ) ;
key_future_time = min ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) >= 0 ) ) ;
if isempty(key_past_time) == 1
    key_past_time = key_future_time ;
end

%%%% find keyframe indexes and fractional contribution
key_past_index = find( INTERPSTACK.time == key_past_time ) ;
key_future_index = find( INTERPSTACK.time == key_future_time ) ;
dist_to_past = abs( key_past_time - t_geol ) ;
dist_to_future = abs( key_future_time - t_geol ) ;
%%%% fractional contribution of each keyframe
if dist_to_past + dist_to_future == 0
    contribution_past = 1 ;
    contribution_future = 0 ;
else
    contribution_past = dist_to_future / ( dist_to_past + dist_to_future ) ;
    contribution_future = dist_to_past / ( dist_to_past + dist_to_future ) ;
end

%%%% intrepolate keyframe CO2 concentrations to generate keyframe fields
%%%% find INTERPSTACK keyframes using model CO2
key_upper_CO2 = min ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) >= 0 ) ) ;
key_lower_CO2 = max ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) <= 0 ) ) ;
%%%% find keyframe indexes and fractional contribution
key_upper_CO2_index = find( INTERPSTACK.CO2 == key_upper_CO2 )  ;
key_lower_CO2_index = find( INTERPSTACK.CO2 == key_lower_CO2 ) ;
dist_to_upper = abs( key_upper_CO2 - CO2ppm ) ;
dist_to_lower = abs( key_lower_CO2 - CO2ppm ) ;

%%%% fractional contribution of each keyframe
if dist_to_upper + dist_to_lower == 0
    contribution_lower = 1 ;
    contribution_upper = 0 ;
else
    contribution_upper = dist_to_lower / ( dist_to_upper + dist_to_lower ) ;
    contribution_lower = dist_to_upper / ( dist_to_upper + dist_to_lower ) ;
end


load('insolation_and_temperature/latitude_inso_99_89.mat')


switch Insolation_mode
    case 1
        % copying longitude band times of latitude-varied insolation lists into an matrice
        if t_geol >=  -524.3 & t_geol <= -514.3
            % Only use the complete Insolation forces 
            ETP = inso_forces(t_geol);
            insolation = inso_forces(t_geol);
        
        else
            insolation = zeros(40,48);
            ETP = zeros(40,1);
        end
    case 2
        if t_geol >=  -524.3 & t_geol <= -514.3
            % Use the high envelope of the insolation forces
            ETP = inso2(t_geol);
            [highenv,] = envelope(latitude_inso_99_89,10,'peak');
            insolation = (ETP - mean(highenv)') .* ones(40,48);
        else
            insolation = zeros(40,48);
            ETP = zeros(40,1);
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create time keyframes using CO2 keyfield contributions   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Runoff
RUNOFF_past = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_past_index); 
RUNOFF_future = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_future_index); 
%%%% Tair
Tair_past = contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_past_index); 
Tair_future =  contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_future_index); 

switch Insolation_mode
    case 1
        % dT/dF K/Wm-2， Huyber et al. 2008
        dT_F_past = 0.97;
        dT_F_future = 0.97;
        factor1 = factor1;
    case 2
        % dT_F is settle as a long-term factor, for example, 1.1 W/m-2
        dT_F = 1.1;
        dT_F_past = dT_F;
        dT_F_future = dT_F;
        factor1 = factor1;
end

delta_temp_past = insolation .* dT_F_past .* factor1;
delta_temp_future = insolation .* dT_F_future .* factor1;

Tair_past = Tair_past + delta_temp_past;
Tair_future = Tair_future + delta_temp_future;

%%%% time kayframes that don't depend on CO2
%%%% Topography
TOPO_past = INTERPSTACK.topo(:,:,key_past_index) ; 
TOPO_future = INTERPSTACK.topo(:,:,key_future_index) ; 

%%%% last keyframe land recorded for plot
land_past = INTERPSTACK.land(:,:,key_past_index) ; 
land_future = INTERPSTACK.land(:,:,key_future_index) ; 

%%%% gridbox area
GRID_AREA_km2 = INTERPSTACK.gridarea ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial silicate weathering   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% West / Maffre weathering approximation

%%%% runoff in mm/yr
% Q_past = RUNOFF_past + delta_temp * 4.167;
% Q_past = RUNOFF_past + delta_temp * 8.141;
Q_past = RUNOFF_past + delta_temp_past * 26.032;
% Q_past = RUNOFF_past;
Q_past(Q_past<0) = 0 ;
% Q_future = RUNOFF_future + delta_temp * 4.167;
% Q_future = RUNOFF_future + delta_temp * 8.141;
Q_future = RUNOFF_future + delta_temp_future * 26.032;
% Q_future = RUNOFF_future;
Q_future(Q_future<0) = 0 ;

%%%% temp in kelvin
T_past = Tair_past + 273 ;
T_future = Tair_future + 273 ;

%%%% Erosion rate from Maffre et al. 2018
model_elevation_past = TOPO_past ;
model_elevation_future = TOPO_future ;
model_elevation_past(isnan(model_elevation_past)==1) = 0 ;
model_elevation_future(isnan(model_elevation_future)==1) = 0 ;

%%%%% make lat and lon grids for erosion calculation
lat_grid = INTERPSTACK.lat'.*ones(40,48) ;
lon_grid = INTERPSTACK.lon.*ones(40,48) ;
%%%%% gradient calculation
[~,~,dFdyNorth_past,dFdxEast_past] = gradientm(lat_grid,lon_grid,model_elevation_past) ;
[~,~,dFdyNorth_future,dFdxEast_future] = gradientm(lat_grid,lon_grid,model_elevation_future) ;
%%%%% topographic slope
tslope_past = ( dFdyNorth_past.^2 + dFdxEast_past.^2 ).^0.5 ;
tslope_future = ( dFdyNorth_future.^2 + dFdxEast_future.^2 ).^0.5 ;
%%%%% pierre erosion calculation, t/m2/yr
k_erosion = 3.3e-3 ; %%%% for 16Gt present day erosion in FOAM
EPSILON_past = k_erosion .* (Q_past.^0.31) .* tslope_past .* max(Tair_past,2) ;
EPSILON_future = k_erosion .* (Q_future.^0.31) .* tslope_future .* max(Tair_future,2) ;
%%%%% check total tonnes of erosion - should be ~16Gt
EPSILON_per_gridbox_past = EPSILON_past .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
EPSILON_per_gridbox_future = EPSILON_future .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
erosion_tot_past = sum(sum(EPSILON_per_gridbox_past)) ;
erosion_tot_future = sum(sum(EPSILON_per_gridbox_future)) ;
erosion_tot = erosion_tot_past*contribution_past + erosion_tot_future*contribution_future ;

%%%% Pierre weathering equation params
Xm = 0.1 ;
K = 6e-5 ; 
kw = 1e-3 ;
Ea = 20 ; 
z = 10 ; 
sigplus1 = 0.9 ; 
T0 = 286 ;
R = 8.31e-3 ;

%%%% equations
R_T_past = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_past) ) ) ;
R_T_future = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_future) ) ) ;
R_Q_past = 1 - exp( -1.*kw .* Q_past ) ;
R_Q_future = 1 - exp( -1.*kw .* Q_future ) ;
R_reg_past = ( (z./EPSILON_past).^sigplus1 ) ./ sigplus1 ;
R_reg_future = ( (z./EPSILON_future).^sigplus1 ) ./ sigplus1 ;

%%%% equation for CW per km2 in each box
CW_per_km2_past = 1e6 .* EPSILON_past .* Xm .* ( 1 - exp( -1.* K .* R_Q_past .* R_T_past .* R_reg_past ) ) ; 
CW_per_km2_future = 1e6 .* EPSILON_future .* Xm .* ( 1 - exp( -1.* K .* R_Q_future .* R_T_future .* R_reg_future ) ) ; 
%%%% CW total
CW_past = CW_per_km2_past .* GRID_AREA_km2 ;
CW_future = CW_per_km2_future .* GRID_AREA_km2 ;

%%%% world CW
CW_past(isnan(CW_past)==1) = 0 ;
CW_future(isnan(CW_future)==1) = 0 ;

% only suits to Cambrian, since no high latitudes in northern hemesphiere
CW_sum_future_low_lat = sum(sum(CW_future(1:32,:)));
CW_sum_future_high_lat = sum(sum(CW_future(33:end,:)));

CW_sum_past_low_lat = sum(sum(CW_past(1:32,:)));
CW_sum_past_high_lat = sum(sum(CW_past(33:end,:)));

CW_tot_low_lat = CW_sum_past_low_lat*contribution_past + CW_sum_future_low_lat*contribution_future;
CW_tot_high_lat = CW_sum_past_high_lat*contribution_past + CW_sum_future_high_lat*contribution_future;

CW_sum_past = sum(sum(CW_past)) ;
CW_sum_future = sum(sum(CW_future)) ;
CW_tot = CW_sum_past*contribution_past + CW_sum_future*contribution_future ;

%%%% carbonate weathering spatial approximation, linear with runoff
k_carb_scale = 200 ; %%%% scaling parameter to recover present day rate
CWcarb_per_km2_past = k_carb_scale * Q_past ; 
CWcarb_per_km2_future = k_carb_scale * Q_future ; 
%%%% CW total
CWcarb_past = CWcarb_per_km2_past .* GRID_AREA_km2 ;
CWcarb_future = CWcarb_per_km2_future .* GRID_AREA_km2 ;


%%%% world CWcarb
CWcarb_past(isnan(CWcarb_past)==1) = 0 ;
CWcarb_future(isnan(CWcarb_future)==1) = 0 ;
CWcarb_sum_past = sum(sum(CWcarb_past)) ;
CWcarb_sum_future = sum(sum(CWcarb_future)) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Grid interpolated variables   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% silicate weathering scale factor by present day rate in cation tonnes
silw_scale = 4.2e8 ; %%%% for k erosion 3.3e-3
%%%% overall spatial weathering
silw_spatial = CW_tot * ( (pars.k_basw + pars.k_granw) / silw_scale) ;
carbw_spatial = ( CWcarb_sum_past*contribution_past + CWcarb_sum_future*contribution_future ) ;

%%%% global average surface temperature
GAST = mean(mean( Tair_past .* pars.rel_contrib ))*contribution_past  +  mean(mean( Tair_future .* pars.rel_contrib ))*contribution_future  ;

delta_temp = delta_temp_past .* contribution_past + delta_temp_future .* contribution_future;

%%%% set assumed ice temperature
Tcrit = -10 ;
%%%% ice line calculations
Tair_past_ice = Tair_past ;
Tair_past_ice(Tair_past_ice >= Tcrit) = 0 ;
Tair_past_ice(Tair_past_ice < Tcrit) = 1 ;
Tair_future_ice = Tair_future ;
Tair_future_ice(Tair_future_ice >= Tcrit) = 0 ;
Tair_future_ice(Tair_future_ice < Tcrit) = 1 ;
%%%% count only continental ice
Tair_past_ice = Tair_past_ice.* land_past ;
Tair_future_ice = Tair_future_ice.* land_future ; 
%%%% sum into lat bands
latbands_past = sum(Tair_past_ice,2) ;
latbands_future = sum(Tair_future_ice,2) ;
latbands_past(latbands_past>0) = 1 ;
latbands_future(latbands_future>0) = 1 ;
%%%% find appropiate lat
latresults_past =  INTERPSTACK.lat .* latbands_past' ;
latresults_future =  INTERPSTACK.lat .* latbands_future' ;
latresults_past(latresults_past == 0) = 90 ;
latresults_future(latresults_future == 0) = 90 ;
%%%% lowest glacial latitude
iceline_past = min(abs(latresults_past)) ;
iceline_future = min(abs(latresults_future)) ;
iceline = iceline_past * contribution_past +  iceline_future * contribution_future ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Global variables   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% effect of temp on VEG %%%% fixed
V_T = 1 - (( (GAST - 25)/25 )^2) ;

%%%% effect of CO2 on VEG
P_atm = CO2atm*1e6 ;
P_half = 183.6 ;
P_min = 10 ;
V_co2 = (P_atm - P_min) / (P_half + P_atm - P_min) ;

%%%% effect of O2 on VEG
V_o2 = 1.5 - 0.5*(O/pars.O0) ; 

%%%% full VEG limitation
V_npp = 2*EVO*V_T*V_o2*V_co2 ;

%%%% COPSE reloaded fire feedback
ignit = min(max(48*mrO2 - 9.08 , 0) , 5 ) ;
firef = pars.kfire/(pars.kfire - 1 + ignit) ;

%%%%% Mass of terrestrial biosphere
VEG = V_npp * firef ;

%%%%%% basalt and granite temp dependency - direct and runoff
Tsurf = GAST + 273 ;
TEMP_gast = Tsurf ;

%%%% COPSE reloaded fbiota
V = VEG ;
f_biota = ( 1 - min( V*W , 1 ) ) * PREPLANT * (RCO2^0.5) + (V*W) ;

%%%% version using gran area and conserving total silw
basw = silw_spatial * ( pars.basfrac * BAS_AREA / ( pars.basfrac * BAS_AREA + ( 1 - pars.basfrac  ) *GRAN_AREA ) ) ;
granw = silw_spatial * ( ( 1 - pars.basfrac  ) *GRAN_AREA / ( pars.basfrac * BAS_AREA + ( 1 - pars.basfrac  ) *GRAN_AREA ) ) ;

%%%% add fbiota
basw = basw * f_biota ;
granw = granw * f_biota ;
carbw = carbw_spatial * f_biota ;

%%%% overall weathering
silw = basw + granw ;
carbw_relative = (carbw/pars.k_carbw) ;

%%%% oxidative weathering 
oxidw = pars.k_oxidw*carbw_relative*(G/pars.G0)*((O/pars.O0)^pars.a);



%%%% pyrite weathering 
pyrw = pars.k_pyrw*carbw_relative*(PYR/pars.PYR0) ;

%%%% gypsum weathering 
gypw = pars.k_gypw*(GYP/pars.GYP0)*carbw_relative ;

%%%%% seafloor weathering, revised following Brady and Gislason but not directly linking to CO2
f_T_sfw = exp(0.0608*(Tsurf-288)) ; 
sfw = pars.k_sfw * f_T_sfw * DEGASS ; %%% assume spreading rate follows degassing here

%%%%%%% Degassing 
ocdeg = pars.k_ocdeg*DEGASS*(G/pars.G0);
ccdeg = pars.k_ccdeg*DEGASS*(C/pars.C0)*Bforcing ;
pyrdeg = pars.k_pyrdeg*(PYR/pars.PYR0)*DEGASS;
gypdeg = pars.k_gypdeg*(GYP/pars.GYP0)*DEGASS;

%%%% gypsum burial
% mgsb = pars.k_mgsb*(S/pars.S0);
mgsb = pars.k_mgsb*(S/pars.S0)*(1/SHORELINE) ;

%%%% carbonate burial
mccb = carbw + silw ;

%epsilon = interp1qr([-1000 -466 -444 -411 -399 0 ]',[0.55 0.55 1.5 1.5 1 1 ]',t_geol) ;
epsilon = interp1qr([-1000 -466 -444 -411 -399 0 ]',[0.5 0.5 1.5 1.5 1 1 ]',t_geol) ;
%%%% COPSE reloaded P weathering
pfrac_silw = 0.8 ;
pfrac_carbw = 0.14 ;
pfrac_oxidw = 0.06 ;
phosw = epsilon * pars.k_phosw*( (pfrac_silw)*( silw/pars.k_silw )  +   (pfrac_carbw)*( carbw/pars.k_carbw ) +  (pfrac_oxidw)*(  oxidw/ pars.k_oxidw )  )  ;

%%%% COPSE reloaded
pland = pars.k_landfrac * VEG * phosw  ;
pland0 = pars.k_landfrac*pars.k_phosw;
psea = phosw - pland ;



%%%% convert total reservoir moles to micromoles/kg concentration    
%Pconc = ( P/pars.P0 ) * 2.2 ;
Pconc_i = ( P_i/pars.P0_i ) * 2.2 ;
Pconc_s = ( P_s/pars.P0_s ) * 2.2 ;
%Nconc = ( N/pars.N0 ) * 30.9 ;
Nconc_i = ( N_i/pars.N0_i ) * 30.9 ;
Nconc_s = ( N_s/pars.N0_s ) * 30.9 ;
%newp = 117 * min(Nconc/16,Pconc) ; 
newp_i = 117 * min(Nconc_i/16,Pconc_i) ; 
newp_s = 117 * min(Nconc_s/16,Pconc_s) ; 

%%%%% carbon burial
%mocb = pars.k_mocb*((newp/pars.newp0)^pars.b) * CB ;
mocb_i = pars.k_mocb_i *((newp_i/pars.newp0_i)^pars.b) * CB ;
mocb_s = pars.k_mocb_s *((newp_s/pars.newp0_s)^pars.b) * CB ;

locb = pars.k_locb*(pland/pland0)*CPLAND  ;



%%%%%% OCEAN ANOXIC FRACTION
k_anox = 10 ; 
k_u = 0.4 ;
%ANOX = 1 / ( 1 + exp( -1 * k_anox * ( k_u * (newp/pars.newp0) - (O/pars.O0) ) ) ) ;
ANOX_i = 1 / ( 1 + exp( -1 * k_anox * ( k_u * (newp_i/pars.newp0_i) - (O/pars.O0) ) ) ) ;
ANOX_s = 1 / ( 1 + exp( -1 * k_anox * ( k_u * (newp_s/pars.newp0_s) - (O/pars.O0) ) ) ) ;

% PYR burial function (COPSE)
%fox= 1/(O/pars.O0) ; 
%%%% mpsb scales with mocb so no extra uplift dependence
%mpsb = pars.k_mpsb*(S/pars.S0)*fox*(mocb/pars.k_mocb) ;
% mpsb_i = pars.k_mpsb_i*(S/pars.S0)*fox*(mocb_i/pars.k_mocb_i) ;
% mpsb_s = pars.k_mpsb_s*(S/pars.S0)*fox*(mocb_s/pars.k_mocb_s) ;
% mpsb_i = pars.k_mpsb_i*(S/pars.S0)*(0.95 + 0.05 * ANOX_i / (1-pars.k_oxfrac_i))*(mocb_i/pars.k_mocb_i) ; % fourth change
% mpsb_s = pars.k_mpsb_s*(S/pars.S0)*(0.95 + 0.05 * ANOX_s / (1-pars.k_oxfrac_s))*(mocb_s/pars.k_mocb_s) ;
mpsb_i = pars.k_mpsb_i*(S/pars.S0)*(0.94 + 0.06 * ANOX_i / (1-pars.k_oxfrac_i))*(mocb_i/pars.k_mocb_i) ; % fourth change
mpsb_s = pars.k_mpsb_s*(S/pars.S0)*(0.94 + 0.06 * ANOX_s / (1-pars.k_oxfrac_s))*(mocb_s/pars.k_mocb_s) ;

%%%%%%% nutrient burial
CNsea = 37.5 ;
monb_i = (mocb_i)/CNsea ;
monb_s = (mocb_s)/CNsea ;

%%% P burial with bioturbation on
% CPbiot = 250 ;
% CPlam = 1000 ;
% % mopb = mocb*( (f_biot/CPbiot) + ( (1-f_biot)/CPlam ) ) ;
% mopb_i = mocb_i*( (f_biot/CPbiot) + ( (1-f_biot)/CPlam ) ) ;
% mopb_s = mocb_s*( (f_biot/CPbiot) + ( (1-f_biot)/CPlam ) ) ;

k_oxic = 250 ;
k_anoxic = 4000 ;
CP_i = k_oxic*k_anoxic / ( (1-ANOX_i)*k_anoxic + ANOX_i*k_oxic ) ;
CP_s = k_oxic*k_anoxic / ( (1-ANOX_s)*k_anoxic + ANOX_s*k_oxic ) ;
mopb_i = mocb_i / CP_i ;
mopb_s = mocb_s / CP_s ;

%capb = pars.k_capb*( mocb/pars.k_mocb ) ;
% capb_i = pars.k_capb_i*( mocb_i/pars.k_mocb_i ) ;
% capb_s = pars.k_capb_s*( mocb_s/pars.k_mocb_s ) ;
% capb_i = pars.k_capb_i*( mocb_i/pars.k_mocb_i )*(O/pars.O0)  ; %third change
% capb_s = pars.k_capb_s*( mocb_s/pars.k_mocb_s )*(O/pars.O0) ;
capb_i = pars.k_capb_i*( mocb_i/pars.k_mocb_i )*((1-ANOX_i) / pars.k_oxfrac_i) ; % third change
capb_s = pars.k_capb_s*( mocb_s/pars.k_mocb_s )*((1-ANOX_s) / pars.k_oxfrac_s) ;

%%%% reloaded
% %fepb = (pars.k_fepb/pars.k_oxfrac)*(1-ANOX)*(P/pars.P0) ;
% fepb_i = (pars.k_fepb_i/pars.k_oxfrac_i)*(1-ANOX_i)*(P_i/pars.P0_i) ;
% fepb_s = (pars.k_fepb_s/pars.k_oxfrac_s)*(1-ANOX_s)*(P_s/pars.P0_s) ;
fepb_i = (pars.k_fepb_i)*(P_i/pars.P0_i)*( 0.50 + 0.50 * (1-ANOX_i) / pars.k_oxfrac_i) ;
fepb_s = (pars.k_fepb_s)*(P_s/pars.P0_s)*( 0.50 + 0.50 * (1-ANOX_s) / pars.k_oxfrac_s) ;

%%%%% nitrogen cycle
%%%% COPSE reloaded
if (N_i/16) < P_i
    nfix_i = pars.k_nfix_i * ( ( ( P_i - (N_i/16)  ) / (  pars.P0_i - (pars.N0_i/16)    ) )^2 ) ;
else
    nfix_i = 0 ;
end

if (N_s/16) < P_s
    nfix_s = pars.k_nfix_s * ( ( ( P_s - (N_s/16)  ) / (  pars.P0_s - (pars.N0_s/16)    ) )^2 ) ;
else
    nfix_s = 0 ;
end

% nfix = nfix_i + nfix_s;

% denit = pars.k_denit * ( 1 + ( ANOX / (1-pars.k_oxfrac) )  ) * (N/pars.N0) ;
% denit_i = pars.k_denit_i * ( 1 + ( ANOX_i / (1-pars.k_oxfrac_i) )  ) * (N_i/pars.N0_i) ;
% denit_s = pars.k_denit_s * ( 1 + ( ANOX_s / (1-pars.k_oxfrac_s) )  ) * (N_s/pars.N0_s) ;
denit_i = pars.k_denit_i * ( ANOX_i / (1-pars.k_oxfrac_i) )   * (N_i/pars.N0_i) ;
denit_s = pars.k_denit_s * ( ANOX_s / (1-pars.k_oxfrac_s) )   * (N_s/pars.N0_s) ;

% denit = denit_i + denit_s;

%%%% reductant input
reductant_input = pars.k_reductant_input * DEGASS ;

%%%%% P and N transfer fluxes (moles/yr) [reactive P (soluble reactive phosphorus) and nitrate]
k_p_s_i = (3.94e10) / pars.P0_s;
p_transfer_shelf_i = P_s * k_p_s_i  ; %present day 3.71E11 moles/yr
%p_transfer_shelf_i = 3.71E1 ; %present day 3.71E11 moles/yr
k_p_i_s = (3.14e10) / pars.P0_i ;
p_transfer_i_shelf =  P_i * k_p_i_s ; %  present day 3.59E11 moles/yr
%p_transfer_i_shelf = 3.59E11;

k_n_s_i = (5.84e12) / pars.N0_s  ; % present day 3.9e14 / 14
n_transfer_shelf_i = N_s * k_n_s_i  ; 

k_n_i_s = (5.84e12) / pars.N0_i ;
%n_transfer_i_shelf = k_n_i_s * N_i ; % present day 4.5e14 /14
n_transfer_i_shelf = N_i * k_n_i_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir calculations  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Interior Ocean Phosphate
dy(1) = p_transfer_shelf_i - p_transfer_i_shelf - mopb_i - capb_i - fepb_i  ;

%%% Shelf Phosphate
dy(22) = psea - p_transfer_shelf_i + p_transfer_i_shelf - mopb_s - capb_s - fepb_s ;

%%% Oxygen
dy(2) = locb + mocb_i + mocb_s - oxidw  - ocdeg  + 2*(mpsb_i+mpsb_s - pyrw  - pyrdeg) - reductant_input ;

%%% Carbon dioxide
dy(3) = -locb - mocb_i- mocb_s + oxidw + ocdeg + ccdeg + carbw - mccb - sfw  + reductant_input ;

%%% Sulphate
dy(4) = gypw + pyrw - mgsb - mpsb_i - mpsb_s + gypdeg + pyrdeg ;

%%%Buried organic C
dy(5) = locb + mocb_s + mocb_i - oxidw - ocdeg ;

%%% Buried carb C 
dy(6) = mccb + sfw - carbw - ccdeg ;

%%% Buried pyrite S
dy(7) = mpsb_s + mpsb_i - pyrw - pyrdeg ;

%%% Assume that buried gypsum S don't vary over time in the Cambrian
dy(8) = 0 ;

%%%% Nitrate ocean interior
dy(11) = n_transfer_shelf_i - n_transfer_i_shelf + nfix_i - denit_i - monb_i;

%%%% Nitrate shelf
dy(23) = - n_transfer_shelf_i + n_transfer_i_shelf + nfix_s - denit_s - monb_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Isotope reservoirs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% d13c and d34s for forwards model
d13c_A = y(16) / y(3) ;
d34s_S = y(17) / y(4) ;

%%%% carbonate fractionation
delta_locb = d13c_A - capdelC_land ; 
delta_mocb = d13c_A - capdelC_marine ; 
delta_mccb = d13c_A ;

%%%%% S isotopes (copse)
delta_mpsb = d34s_S - capdelS ;

%%% deltaORG_C*ORG_C 
dy(12) =  locb*(  delta_locb ) + (mocb_i + mocb_s)*( delta_mocb )  -   oxidw*delta_G  -   ocdeg*delta_G ;

%%%%%% deltaCARB_C*CARB_C 
dy(13) =  mccb*delta_mccb + sfw*delta_mccb  -  carbw*delta_C  - ccdeg*delta_C  ;

%%% deltaPYR_S*PYR_S (young)
dy(14) =  (mpsb_i + mpsb_s)*( delta_mpsb )  - pyrw*delta_PYR  - pyrdeg*delta_PYR  ;

%%% deltaGYP_S*GYP_S (young)
dy(15) =  mgsb*d34s_S    - gypw*delta_GYP  - gypdeg*delta_GYP  ;

%%% delta_A * A
dy(16) = -locb*(  delta_locb ) -(mocb_s + mocb_i)*( delta_mocb ) + oxidw*delta_G + ocdeg*delta_G + ccdeg*delta_C + carbw*delta_C - mccb*delta_mccb - sfw*delta_mccb + reductant_input*-5 ;

%%% delta_S * S
dy(17) = gypw*delta_GYP + pyrw*delta_PYR -mgsb*d34s_S - (mpsb_i + mpsb_s)*( delta_mpsb ) + gypdeg*delta_GYP + pyrdeg*delta_PYR ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Strontium system   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% fluxes
Sr_granw = pars.k_Sr_granw *( granw / pars.k_granw ) ;
Sr_basw = pars.k_Sr_basw *( basw / pars.k_basw ) ;
Sr_sedw = pars.k_Sr_sedw *( carbw / pars.k_carbw ) * (SSr/pars.SSr0) ;
Sr_mantle = pars.k_Sr_mantle * DEGASS ;
Sr_sfw = pars.k_Sr_sfw * (sfw/pars.k_sfw) * ( OSr/pars.OSr0 ) ;
Sr_metam = pars.k_Sr_metam * DEGASS * (SSr/pars.SSr0) ;
Sr_sedb = pars.k_Sr_sedb * ( mccb/pars.k_mccb ) * ( OSr/pars.OSr0 ) ;

%%%% fractionation calculations
delta_OSr = y(19) / y(18) ;
delta_SSr = y(21) / y(20) ;

%%%% original frac
RbSr_bas = 0.1 ;
RbSr_gran = 0.26 ;
RbSr_mantle = 0.066 ;
RbSr_carbonate = 0.5 ;

%%%% frac calcs
dSr0 = 0.69898 ;
tforwards = 4.5e9 + t ;
lambda = 1.4e-11 ;
dSr_bas = dSr0 + RbSr_bas*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_gran = dSr0 + RbSr_gran*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_mantle = dSr0 + RbSr_mantle*( 1 - exp(-1*lambda*tforwards) ) ;

%%%% Ocean [Sr]
dy(18) = Sr_granw + Sr_basw + Sr_sedw + Sr_mantle - Sr_sedb - Sr_sfw ;

%%%% Ocean [Sr]*87/86Sr
dy(19) = Sr_granw*dSr_gran + Sr_basw*dSr_bas + Sr_sedw*delta_SSr + Sr_mantle*dSr_mantle - Sr_sedb*delta_OSr - Sr_sfw*delta_OSr ;

%%%% Sediment [Sr]
dy(20) = Sr_sedb - Sr_sedw - Sr_metam ;

%%%% Sediment [Sr]*87/86Sr
dy(21) = Sr_sedb*delta_OSr - Sr_sedw*delta_SSr - Sr_metam*delta_SSr + SSr*lambda*RbSr_carbonate*exp(lambda*tforwards)  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Mass conservation check   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_C = A + G + C ;
res_S = S + PYR + GYP ;
iso_res_C = A*d13c_A + G*delta_G + C*delta_C ;
iso_res_S = S*d34s_S + PYR*delta_PYR + GYP*delta_GYP ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record full states for single run   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensanal == 0
    workingstate.iso_res_C(stepnumber,1) = iso_res_C ;
    workingstate.iso_res_S(stepnumber,1) = iso_res_S ;
    workingstate.res_C(stepnumber,1) = res_C ;
    workingstate.res_S(stepnumber,1) = res_S ;
    workingstate.time(stepnumber,1) = t;
    workingstate.temperature(stepnumber,1) = TEMP_gast ;
    workingstate.tempC(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.P_i(stepnumber,1) = P_i ;
    workingstate.P_s(stepnumber,1) = P_s ;
    workingstate.O(stepnumber,1) = O ;
    workingstate.A(stepnumber,1) = A ;
    workingstate.S(stepnumber,1) = S ;
    workingstate.G(stepnumber,1) = G ;
    workingstate.C(stepnumber,1) = C ;
    workingstate.PYR(stepnumber,1) = PYR ;
    workingstate.GYP(stepnumber,1) = GYP ;
    workingstate.N_i(stepnumber,1) = N_i ;
    workingstate.N_s(stepnumber,1) = N_s ;
    workingstate.OSr(stepnumber,1) = OSr ;
    workingstate.SSr(stepnumber,1) = SSr ;
    %%%%%%% print isotope information
    workingstate.d13c_A(stepnumber,1) = d13c_A ;
    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_G(stepnumber,1) = delta_G ;
    workingstate.delta_C(stepnumber,1) = delta_C ;
    workingstate.delta_PYR(stepnumber,1) = delta_PYR ;
    workingstate.delta_GYP(stepnumber,1) = delta_GYP ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    %%%%%%% print forcings
    workingstate.DEGASS(stepnumber,1) = DEGASS ;
    workingstate.W(stepnumber,1) = W ;
    workingstate.EVO(stepnumber,1) = EVO ;
    workingstate.CPLAND(stepnumber,1) = CPLAND ;
    workingstate.Bforcing(stepnumber,1) = Bforcing ;
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA ;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA ;
    %%%%%%%% print variables
    workingstate.RCO2(stepnumber,1) = RCO2 ;
    workingstate.RO2(stepnumber,1) = RO2 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.VEG(stepnumber,1) = VEG ;
    workingstate.ANOX_s(stepnumber, 1) = ANOX_s ; 
    workingstate.ANOX_i(stepnumber, 1) = ANOX_i ;
    workingstate.iceline(stepnumber,1) = iceline ;
    %%%%%%%% print fluxes
    workingstate.mocb_i(stepnumber,1) = mocb_i ;
    workingstate.mocb_s(stepnumber,1) = mocb_s ;
    workingstate.locb(stepnumber,1) = locb ;
    workingstate.mccb(stepnumber,1) = mccb ;
    workingstate.mpsb_i(stepnumber,1) = mpsb_i ;
    workingstate.mpsb_s(stepnumber,1) = mpsb_s ;
    workingstate.mgsb(stepnumber,1) = mgsb ;
    workingstate.silw(stepnumber,1) = silw ;
    workingstate.carbw(stepnumber,1) = carbw ;
    workingstate.oxidw(stepnumber,1) = oxidw ;
    workingstate.basw(stepnumber,1) = basw ;
    workingstate.granw(stepnumber,1) = granw ;
    
    workingstate.phosw(stepnumber,1) = phosw ;
    workingstate.psea(stepnumber,1) = psea ;
    workingstate.fepb_i(stepnumber,1) = fepb_i ;
    workingstate.fepb_s(stepnumber,1) = fepb_s ;
    workingstate.capb_i(stepnumber,1) = capb_i ;
    workingstate.capb_s(stepnumber,1) = capb_s ;
    workingstate.mopb_i(stepnumber,1) = mopb_i ;
    workingstate.mopb_s(stepnumber,1) = mopb_s ;
    
    workingstate.nfix_i(stepnumber,1) = nfix_i ;
    workingstate.nfix_s(stepnumber,1) = nfix_s ;
    workingstate.denit_i(stepnumber,1) = denit_i ;
    workingstate.denit_s(stepnumber,1) = denit_s ;
    workingstate.VEG(stepnumber,1) = VEG ;
    workingstate.pyrw(stepnumber,1) = pyrw ;
    workingstate.gypw(stepnumber,1) = gypw ;
    workingstate.ocdeg(stepnumber,1) = ocdeg ;
    workingstate.ccdeg(stepnumber,1) = ccdeg ;
    workingstate.pyrdeg(stepnumber,1) = pyrdeg ;
    workingstate.gypdeg(stepnumber,1) = gypdeg ;
    workingstate.sfw(stepnumber,1) = sfw ;
    workingstate.Sr_granw(stepnumber,1) = Sr_granw ;
    workingstate.Sr_basw(stepnumber,1) = Sr_basw ;
    workingstate.Sr_sedw(stepnumber,1) = Sr_sedw ;
    workingstate.Sr_mantle(stepnumber,1) = Sr_mantle ;
    workingstate.dSSr(stepnumber,1) = dSSr ;
    workingstate.relativenewp_i(stepnumber,1) = newp_i/pars.newp0_i ;
    workingstate.relativenewp_s(stepnumber,1) = newp_s/pars.newp0_s ;
    workingstate.erosion_tot(stepnumber,1) = erosion_tot ;
    
    %%% insolation matrix save
    workingstate.ETP(stepnumber,:) = ETP' ;
    workingstate.delta_temp(stepnumber,:) = delta_temp(:,1)' ;

    %%%save the low/high latitude and total continental weathering flux
    workingstate.CW_tot(stepnumber,1) = CW_tot;
    workingstate.CW_tot_low_lat(stepnumber,1) = CW_tot_low_lat;
    workingstate.CW_tot_high_lat(stepnumber,1) = CW_tot_high_lat;
    %%%%%%%% print time
    if pars.runcontrol == -2
        workingstate.time_myr(stepnumber,1) = t.*1e-6 ;
    else
        workingstate.time_myr(stepnumber,1) = t_geol ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record gridstates for single run   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 0
    %%%% print a gridstate when each keytime threshold is crossed, or at model end
    next_stamp = pars.next_gridstamp ;
    if pars.finishgrid == 0
        if t_geol > next_stamp || t_geol == 0

            %%%% write gridstates
            gridstate.time_myr(pars.gridstamp_number,1) = next_stamp ;
            gridstate.land(:,:,pars.gridstamp_number) = land_past ;
            gridstate.Q(:,:,pars.gridstamp_number) = Q_past ;
            gridstate.Tair(:,:,pars.gridstamp_number) = Tair_past ;
            gridstate.TOPO(:,:,pars.gridstamp_number) = TOPO_past ;
            gridstate.CW(:,:,pars.gridstamp_number) = CW_per_km2_past ; %%% t/km2/yr
            gridstate.CWcarb(:,:,pars.gridstamp_number) = CWcarb_past ;
            gridstate.EPSILON(:,:,pars.gridstamp_number) = EPSILON_past * 1e6 ; %%% t/km2/yr

            %%%% set next boundary
            if t_geol < 0
                pars.gridstamp_number = pars.gridstamp_number + 1 ;
                pars.next_gridstamp = pars.runstamps(pars.gridstamp_number) ;
            else
                pars.finishgrid = 1 ;
            end

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Record plotting states only in sensanal   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 1
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA;
    workingstate.DEGASS(stepnumber,1) = DEGASS;
    workingstate.phosw(stepnumber,1) = phosw ;


    workingstate.ANOX_i(stepnumber,1) = ANOX_i ;
    workingstate.ANOX_s(stepnumber,1) = ANOX_s ;

    workingstate.mocb_i(stepnumber,1) = mocb_i;
    workingstate.mocb_s(stepnumber,1) = mocb_s;

    workingstate.mpsb_s(stepnumber,1) = mpsb_s;
    workingstate.mpsb_i(stepnumber,1) = mpsb_i;

    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    workingstate.SmM(stepnumber,1) = 28*S/pars.S0 ;
    workingstate.CO2ppm(stepnumber,1) = RCO2*280 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.iceline(stepnumber,1) = iceline ;

    workingstate.T_gast(stepnumber,1) = TEMP_gast - 273 ;

    
    workingstate.P_i(stepnumber,1) = P_i/pars.P0_i ;
    workingstate.P_s(stepnumber,1) = P_s/pars.P0_s ;
    workingstate.N_s(stepnumber,1) = N_s/pars.N0_s ;
    workingstate.N_i(stepnumber,1) = N_i/pars.N0_i ;
    workingstate.time_myr(stepnumber,1) = t_geol ;
    workingstate.time(stepnumber,1) = t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Final actions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% output timestep if specified
if sensanal == 0
    if pars.telltime ==1
        if mod(stepnumber,pars.display_resolution) == 0 
            %%%% print model state to screen
            fprintf('Model step: %d \t', stepnumber); fprintf('time: %d \t', t_geol) ; fprintf('next keyframe: %d \n', next_stamp)
        end
    end
end


%%%% record current model step
stepnumber = stepnumber + 1 ;


%%%% option to bail out if model is running aground
if stepnumber > pars.bailnumber
   terminateExecution
end


%%%% end function
end



