% Climate_Sys.m
function [ dxdt , y ] = Climate_Sys( year , x, u, c );
% the ODE's of the EnRoads climate module

%% CONSTANTS via elements of the structure ... c.xyz ---------------------

%% STATES ----------------------------------------------------------------

Q  = x(1:5);  % heat in surface and ocean layers [ Watt year / sq.m ] 
C  = x(6);    %  C  in atmosphere  [ Gton ]
N  = x(7);    % N2O in atmosphere  [ Mton ]
M  = x(8);    % CH4 in atmosphere  [ Mton ]
% PFCs in atm = x(9);   
% SF6 in atm  = x(10);  
% HFCs in atm = x(11);  % stock of each HFC in atmosphere 
% ODS in atm = x(12);   % aggregate stock of CFCs and HCFCs etc regulated by Montreal Protocol 

%% ALGEBRAIC EQUATIONS ---------------------------------------------------

% total GHG radiative forcing -----------------------------------------

% IPCC AR6. WG1. Table 7.SM.1.
% https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf

co2_ppm = C * c.ppm_co2_per_GtonC    / c.unit_ppm;
ch4_ppb = M * c.ppb_ch4_per_mton_ch4 / c.unit_ppb; 
n2o_ppb = N * c.ppb_n2o_per_mton_n2o / c.unit_ppb; 

if (co2_ppm < c.co2_ref)
  alpha_co2 = c.d1;
end
if (co2_ppm > c.co2_ref )
  alpha_co2 = c.d1 + c.a1*(co2_ppm-c.co2_ref)^2 + c.b1*(co2_ppm-c.co2_ref);
end
if (co2_ppm > c.co2_alpha_max)
  alpha_co2 = c.d1 - c.b1^2/(4*c.a1);
end

alpha_n2o = c.c1*sqrt( n2o_ppb );

rf_co2 = ( alpha_co2 + alpha_n2o ) * log( C / c.c_preindustrial);

rf_ch4 = ( c.a3*sqrt(ch4_ppb) + c.b3*sqrt(n2o_ppb) + c.d3 )*( sqrt(ch4_ppb) - sqrt(c.ch4_ref) );

rf_n2o = ( c.a2*sqrt(co2_ppm) + c.b2*sqrt(n2o_ppb) + c.c2*sqrt(ch4_ppb) + c.d2 )*( sqrt(n2o_ppb) - sqrt(c.n2o_ref) );

rf_f_gasses = 0; % ... just for now ...
rf_ods      = 0; % ... just for now ...

well_mixed_ghg_forcing = rf_co2 + rf_ch4 + rf_n2o + rf_f_gasses + rf_ods;

other_forcing = -1.0 + (year - 1990) * ( -0.3 - (-1) ) / (2100 - 1990 );

hydrogen_released_per_year_in_Tg = 0.0; % leave as a constant ... just for now ...
h2_effect_on_rf = hydrogen_released_per_year_in_Tg * ...
( c.ch4_erf_per_h2_flux + c.o3_erf_per_h2_flux + c.strat_h2o_erf_per_h2_flux + c.aerosol_rf_per_h2_flux);

total_radiative_forcing = well_mixed_ghg_forcing + other_forcing + h2_effect_on_rf;

% temperature at all layers -------------------------------------------

temperature = Q ./ c.heat_capacity;

% feedback cooling ----------------------------------------------------

x2_co2_forcing = (alpha_co2 + alpha_n2o) * log(2);

climate_feedback_param = x2_co2_forcing / c.climate_sensitivity_to_2x_co2;

feedback_cooling = temperature(1) * climate_feedback_param;

heat_transfer_ocean_layers = ( temperature(1:4) - temperature(2:5) ) .* c.heat_transfer_coeff ./ ((c.layer_depth(1:4) + c.layer_depth(2:5))/2);


%% DIFFERENTIAL EQUATIONS  -----------------------------------------------

dQdt = [ total_radiative_forcing - feedback_cooling - heat_transfer_ocean_layers(1); % surface
         heat_transfer_ocean_layers(1:4) - [ heat_transfer_ocean_layers(2:4) ; 0 ] ]; % ocean

dCdt =   +0.005*C + 0.000*randn;   % ... just for now ...
dNdt =   +0.004*N + 0.000*randn;   % ... just for now ...
dMdt =   +0.003*M + 0.000*randn;   % ... just for now ...

dxdt = [ dQdt ; 
         dCdt ;
         dNdt ;
         dMdt ];


%% OUTPUTS ---------------------------------------------------------------

y = [ temperature ]; % temperatures in atmosphere and ocean from 1750 deg C

% ==========================================================================
