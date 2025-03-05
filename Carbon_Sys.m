% Carbon_Sys.m
function [ dxdt , y ] = Carbon_Sys( year , x , u , c )
% the ODE's of the EnRoads carbon module (and other GHG's)

%% CONSTANTS via elements of the structure ... c.xyz

%% STATES ----------------------------------------------------------------

C = x(1:6);  %  C  in atmosphere in atmosphere and ocean Gton
N = x(7);    % N2O in atmosphere  Mton
M = x(8);    % CH4 in atmosphere  Mton
% PFCs in atm = x(9);   
% SF6 in atm  = x(10);  
% HFCs in atm = x(11);  % stock of each HFC in atmosphere 
% ODS in atm = x(12);   % aggregate stock of CFCs and HCFCs etc regulated by Montreal Protocol 


%% INPUTS -------------------------------------------------------------

% temperature change from 1850 in atmosphere and oceans deg C
temperature = u - c.mean_temperature_change_from_v_1750_to_v_1850; 


%% ALGEBRAIC RELATIONS ---------------------------------------------------

% Carbon ................................

ch4_uptake = M / c.time_const_for_ch4; %  [ Mtons/Year ]

% C from CH4 oxidation : Flux of C into the atmosphere from the oxidation of CH4, the mode of removal of CH4 from atmosphere.
c_from_ch4_oxidation = ch4_uptake / c.ch4_per_c / c.mton_c_per_gton_c;

% Carbon from emissions, decay, ch4, CDR, ocean, biomass 
%{
% C loss from storage from non land CDR (Carbon Dioxide Removal) 
% Total flow of carbon from geological or other nonLULUCF storage. (Land Use, Land Use Change, Forrestry)
c_loss_from_storage_from_non_land_cdr = total_loss_of_c_from_storage_from_nonlulucf_cdr + total_c_leakage_from_ccs_(carbon_capture_and_storage);

flux_c_biomass_to_atmosphere u(Land use) = lulucf_biomass_emissions_u(Land use) + c_from_biomass_respiration_u[Land use]); 

flux_c_soil_to_atmosphere = sum( flux_c_soil_to_atmosphere_u(Land use) );

flux_c_to_atm_from_structure_decay = sum( flux_c_to_atm_from_industrial_wood_decay_u(Land use) );

global_c_energy_and_industry_emissions = ( global_co2_energy_and_industry_emissions - indirect_co2_in_ch4_accounting ) / co2_per_C;

% Excludes all CO2 emissions from land including those from bioenergy.
... global_co2_energy_and_industry_emissions = IF THEN ELSE( SWITCH SCC = 1, CO2 emissions excluding all LULUCF + 
... 			CO2 emissions excluding all LULUCF  =
... 			SUM( CO2 emissions by emissions source excluding bioenergy[GHG sources] )

% Indirect CO2 in CH4 accounting ..........

% Indirect CO2 emissions included in accounting data that occurs due to oxidation of CH4;
% in this model the indirect emissions are explicit from the methane cycle,
% and therefore are deducted from the CO2 accounting to correct the data.
% Also includes emissions for which CH4 represents recently-extracted biomass
% (e.g. enteric fermentation). 

indirect_co2_in_ch4_accounting = ( ch4_emissions_by_emissions_source(energy_production)
			 + 
			 ch4_emissions_by_emissions_source(End_use_capital)
			 + 
			 ch4_emissions_by_emissions_source(Waste) ) / ch4_per_C_/ mton_c_per_gton_c * co2_per_c * fraction_of_anthro_ch4_emissions_included_in_indirect_co2 accounts;

%C removal by non land CDR : Sequestration via nonLULUCF CDR.  [  GtonsC / Year ]
c_removal_by_non_land_cdr = sum( c_removal_by_nonAF_CDR(Sequestration_type) );

% C removal by nonAF CDR : 
% Annual removal, in Gton carbon for each of the sequestration types using either
% a simple time delay or specific detail structure, with limit based on C in atmosphere.

c_removal_by_nonaf_cdr(nonaf_cdr)  = ...
min( c_removal_rate_for_nonaf_cdr(sequestration_type) , C(1) / length( nonaf_cdr ) / minimum_time_for_c_cycle_changes )

% Flux C Atmosphere to Biomass : Carbon flux from atmosphere to biosphere (from primary production),
% including feedbacks from CO2 fertilization. Calculated as sum of the fluxes for all land types.
% [ GtonsC/Year ]

flux_c_atmosphere_to_biomass = sum( flux_c_atmosphere_to_biomass_u(Land use) );
%}


% Effect of Temperature on DIC pCO2 : (dimensionless)
% The fractional reduction in the solubility of CO2 in ocean falls with rising temperatures.
% We assume a linear relationship, likely a good approximation over the typical range for warming until 2100.
effect_of_temp_on_dic_pco2 = 1.0 - c.sensitivity_of_pco2_dic_to_temperature * temperature(1);

% Equil C in Mixed Layer [ GtonsC ]
% Equilibrium carbon content of mixed layer.
% Determined by the Revelle buffering factor, and by temperature.
% For simplicity, we assume a linear impact of warming on the equilibrium solubility of CO2 in the ocean.
% The user controls the strength of that effect.

% Buffer Factor : Buffer factor for atmosphere/mixed ocean carbon equilibration.  (dimensionless)
buffer_factor = c.ref_buffer_factor * ( C(2) / c.preind_c_in_mixed_layer) ^ c.buff_c_coeff;

equil_c_in_mixed_layer = c.preind_c_in_mixed_layer * effect_of_temp_on_dic_pco2 * ( C(1) / c.c_preindustrial )^(1.0/buffer_factor);

% C Flux Atm to Ocean : Carbon flux from atmosphere to mixed ocean layer.  [ GtonsC/Year ]
flux_c_atm_to_ocean = ( ( equil_c_in_mixed_layer - C(2) ) / c.mixing_time );

% Carbon in the ocean .....................................................

% C in layer_i per meter [ GtonsC / meter ]
c_in_layer_per_meter = C(2:6) ./ c.layer_depth;

% Diffusion Flux : Diffusion flux between ocean layers.  [ GtonsC/Year ]
diffusion_flux = ( c_in_layer_per_meter(1:4) - c_in_layer_per_meter(2:5) ) * c.eddy_diff_coeff ./ ((c.layer_depth(1:4) + c.layer_depth(2:5))/2);

% N2O ...................................

%{
% Total N2O emissions : Name change for the variable that was previously "Global N2O anthro emissions".  [ Mton/Year ]
total_n2o_emissions = sum( n2o_emissions_by_emissions_source );

% Global N2O anthro emissions : [ Mton/Year ]
% Until 2100, anthropogenic N2O emissions follow the normal path set by assumptions and policy scenarios.
% After 2100, they are assumed to change at the same rate as CO2 emissions.
global_n2o_anthro_emissions = total_n2o_emissions;

%Global N2O Emissions : Sum of anthropogenic and natural N2O emissions.  [ Mton/Year ]
global_n2o_emissions = global_n2o_anthro_emissions + c.natural_n2o_emissions;
%}

n2o_uptake = N / c.time_const_for_n2o; %  [ Mton / yr ]

% CH4 ...................................

%{
% Sum of the anthropogenic CH4 as determined from ratio of RS and the additional CH4 from the leakage of natural gas.
total_ch4_emissions = sum( ch4_emissions_by_emissions_source );  % [ Mton / yr ]
global_ch4_anthro_emissions = total_ch4_emissions;
 
% Natural CH4 Emissions : Flux of methane from anaerobic respiration in the biosphere, [ Mtons CH4/year ]
natural_ch4_emissions = ( flux_biosphere_c_to_ch4_natural + flux_ocean_c_to_ch4 ) * ch4_per_c * mton_c_per_gton_c;

% Total Flux Biomass to CH4 natural [ GtonsC/year ]
total_flux_biomass_to_ch4_natural = sum( flux_biomass_c_to_ch4_natural_u )

% Flux Biosphere C to CH4 natural : Carbon flux from biosphere as methane, in GtC/year, arising from anaerobic respiration.
% [ GtonsC/year ]
flux_biosphere_c_to_ch4_natural = total_flux_biomass_to_ch4_natural + flux_soil_c_to_ch4;
%}


%% DIFFERENTIAL EQUATIONS  -----------------------------------------------

% Carbon in Atmosphere and Ocean

dCdt =  [ c_from_ch4_oxidation ...  % atmosphere 
        - flux_c_atm_to_ocean ;  
         % + c_loss_from_storage_from_non_land_cdr
         % + flux_c_biomass_to_atmosphere
         % + flux_c_soil_to_atmosphere
         % + flux_c_to_atm_from_structure_decay
         % + global_c_energy_and_industry_emissions
         % - c_removal_by_non_land_cdr
         % - flux_c_atmosphere_to_biomass ];
          flux_c_atm_to_ocean  ...   %  top   layer  of ocean
          - c.flux_ocean_c_to_ch4 ... 
          - diffusion_flux(1) ;
          diffusion_flux(1:3)  ...   % lower  layers of oceean 
          - diffusion_flux(2:4) ;
          diffusion_flux(4) ];       % bottom layer  of oceean 


% Nitrous Oxide ... N2O in Atmosphere ... [ Mtons ]
dNdt = - n2o_uptake;
     % + global_n2o_emissions 

% Methane ... CH4 in Atmosphere ... [ Mtons ]
dMdt = - ch4_uptake;
     % + global_ch4_anthropomorphic_emissions 
     % + natural_ch4_emissions;

dxdt = [ dCdt ;
         dNdt ;
         dMdt ];

%% OUTPUTS ---------------------------------------------------------
y = 0 ; %[ temperature ];  

% ===================================================================================================
