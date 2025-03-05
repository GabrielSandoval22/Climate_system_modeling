% EnRoads_Constants.m

%% universal constants -------------------------------------------

c.c100_percent = 100;       % [ percent ]
 
c.one_year = 1.0;           % [ year ]

c.unit_ppb = 1.0;           % [ ppb ]
c.unit_ppm = 1.0;           % [ ppm ]

c.hours_per_year = 24*365.25;                % [ hours/year ]
c.seconds_per_year = c.hours_per_year*60*60; % [ sec / year ] 

c.ppt_per_ppb = 1000;       % [ ppt/ppb ]
c.ppt_per_mol = 5.68e-09;   % [ ppt/mole ]

c.g_per_ton = 1e6;          % [ g/ton ]
c.kg_per_tg = 1e+09;        % [ kg / tg ]
c.ton_per_mton = 1e6;       % [ ton / Mton ]
c.mton_c_per_gton_c = 1000; % [ MtonsC/GtonsC ]

c.watt_s_per_J = 1.0;       % definition of a watt [ Watt second per Joule ]
c.MJ_per_EJ = 1e+12;        % [ MJ/EJ ]
c.GJ_per_kWh = 0.0036;      % [ GJ/kWh ]

c.density = 1000;           % density of water [ kg / meter^3 ] 
c.mass_heat_cap = 4186 ;    % heat capacity per unit mass of water [ J/kg/degC ]
% atmosphere and upper ocean heat capacities per unit area
% based on 70% 100m ocean layer and 30% 8.4m equiv land layer [ Watt.year/meter^2 ]
% Heat Capacity per unit mass   (kg)   ... water ...    4 kJ/kg/degC = 4000 Watt sec / kg / deg C  
% Heat Capacity per unit volume (cu.m) ... water ... 4000 kJ/cu.m/degC = 4000000 Wat sec / cu.m / deg C

c.present_year = 2025; 
% Used to shift from data-driven to policy-driven structure.
% There are many instances of policy start time constants with min and default value
% that must be changed manualy when updating at the end of each year. Search for this value in a text view.

%% Earth model constants -------------------------------------------

c.earth_area = 5.1e14 ;        % surface area of the Earth [ meter^2 ]
c.land_area_fraction = 0.292 ; % fraction of earth that is land [ meter^2 ] 

% effective land area heat capacity, expressed as equivalent water layer thickness
c.land_thickness = 8.4; % [ meter ]

% ocean layer thicknesses 
c.layer_depth = [ 100 ; 300 ; 300 ; 1300 ; 1800 ]; % [ meter ]

c.initial_temp_change_from_1750 = [ 0.7 ; 0.34 ; 0.17 ; 0.0156 ; 5.9e-04 ]; % [ deg C ]

c.mean_temperature_change_from_v_1750_to_v_1850 = 0.05; % [ degC ]

% volume of water (and water equivalent) in each layer including surface layer (and atmosphere)
c.layer_volume = c.earth_area * [ c.land_area_fraction * c.land_thickness + (1 - c.land_area_fraction) * c.layer_depth(1)  ; 
                  ( (1 - c.land_area_fraction) * c.layer_depth(2:5) ) ] ; 

c.volumetric_heat_capacity = c.mass_heat_cap * c.watt_s_per_J / c.seconds_per_year * c.density; % [ Watt year / deg C / meter^3 ]

% Specific Heat Capacity: Heat Capacity per unit mass   (kg)   ... water ...    4 kJ/kg/degC
% Specific Heat Capacity: Heat Capacity per unit volume (cu.m) ... water ... 4000 kJ/cu.m/degC
% Rd  : heat capacity of deep ocean
% Rs  : heat capacity of surface water and atmosphere
% tau : heat transfer time constant 

% heat capacity of each layer, including surface layer (and atmosphere)
c.heat_capacity = c.layer_volume .* c.volumetric_heat_capacity / c.earth_area; % [ Watt year / deg C / meter^2 ]

% initial heat in each layer, including surface layer (and atmosphere)
c.heat_initial = c.initial_temp_change_from_1750 .* c.heat_capacity; % [ Watt year / meter^2 ]

c.heat_transfer_rate   = 1.23;      % watt/(meter*meter)/DegreesC
c.heat_diffusion_covar = 1;         % (dimensionless)
%Eddy diff mean : Mean rate of vertical transport at which carbon is mixed in the ocean due to eddy motion. [ meter*meter/year ]
c.eddy_diff_mean = 4400;
c.ocean_heat_mixing_time = 1.0;       % year

% Atmosphere - mixed ocean layer mixing time.
c.mixing_time = 1.0;   % year


% Eddy diffusian coefficient index :
% Index of coefficient for rate at which carbon is mixed in the ocean
% due to eddy motion, where 1 is equivalent to the expected value of 4400 meter/meter/year.
c.eddy_diff_coeff_index = 1; % (dimensionless)
c.eddy_diff_coeff = c.eddy_diff_coeff_index * c.eddy_diff_mean; % [ meter*meter/year ]

% Heat transfer coefficient :
% The vertical transport at which heat is mixed in the ocean due to eddy motion.
% [ watt/(meter*meter)/(DegreesC/meter) ]
c.heat_transfer_coeff = c.heat_transfer_rate * (c.layer_depth(1:4) + c.layer_depth(2:5)) / 2 * c.heat_diffusion_covar * ( ( c.eddy_diff_coeff / c.eddy_diff_mean) + ( 1 - c.heat_diffusion_covar ) );


%% Climate model constants -------------------------------------------
% lambda : Climate Feedback Parameter ... 1 < lambda < 6 ... [ W / sq.m / deg C ] 

c.climate_sensitivity_to_2x_co2 = 3; % deg C
%  Equilibrium temperature change in response to a 2xCO2 equivalent change in radiative forcing.
%  According to AR6, high confidence that value likely between 2.5 and 4,
%  with high confidence that extremely unlikely less than 2 and medium confidence
%  that extremely unlikely greater than 5. Changed from AR5, which had 1.5-4.5 outer limits 1-6.

%  Data from ... 
%  M. Meinshausen, S. Smith et al. 
%  "The RCP GHG concentrations and their extension from 1765 to 2300",
%  Climatic Change. 109(213) (2011)
%  DOI 10.1007/s10584-011-0156-z,
%  https://link.springer.com/article/10.1007/s10584-011-0156-z


%% GHG constants ---------------------------------------------------

c.co2_per_c = 3.66; % [ GtonsCO2 / GtonsC
c.ch4_per_c = 1.33; % GtonsCO2 / GtonsC

% Radiant Heat Forcing Coeffcients from IPCC AR6
% IPCC AR6. WG1. Table 7.SM.1.
% https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf
% RF formulae for CO2, CH4 and N2O.

% CO2 ............................................................
c.a1      =   -2.4785e-7;       % W / m^2 / ppm^2
c.b1      =    7.5906e-4;       % W / m^2 / ppm
c.c1      =   -2.1492e-3;       % W / m^2 / ppb^(0.5)
c.d1      =    5.2488;          % W / m^2
c.C0      =  353.98;            % initial carbon concentration in atmosphere from 1990 Mauna Loa, 2020 [ ppm ] 
c.co2_ref =  277.15;            % CO2 reference concentration [ ppm ] 
% Restated per CDIAC, http://cdiac.ornl.gov/pns/convert.html , that reports:
% 1 ppm by volume of atmosphere CO2 = 2.13 Gt C (Uses atmospheric mass (Ma) = 5.137 x 10^18 kg)
c.ppm_co2_per_GtonC = 0.4695;   % ppm / GtonC
% Preindustrial C [ GtonsC ]
% Calculated from preindustrial concentration of 277 ppm divided by 0.4695 ppm CO2 per GtonsC.
% Atmospheric CO2 record based on ice core data before 1958 (Ethridge et. al., 1996; MacFarling Meure et al., 2006)
c.c_preindustrial = 590;        % preindustrial Carbon concentration (not CO2) ppm 

c.co2_alpha_max = c.C0 - c.b1 / (2*c.a1);  % CO2 threshold in EnRoads = 1808.4 ppm

% N2O ............................................................
% Adjusts total RF from CH4 and N2O to be less than the sum of RF from each individually to account for interactions between both gases.
c.a2      =   -3.4197e-4;       % W / m^2 / ppm
c.b2      =    2.5455e-4;       % W / m^2 / ppb
c.c2      =   -2.4357e-4;       % W / m^2 / ppb
c.d2      =    0.12173;         % W / m^2 / ppb^(0.5)
c.N0      =  308.725;           % initial N2O concentration. based on 1990 levels from GISS, 2020. [ ppb ]
c.n2o_ref =  273.87;            % N2O reference concentration [ ppb ]
c.natural_n2o_emissions = 11.2; % Mton / year
c.n2o_N_molar_mass = 28;        % g / mole  
c.ppb_n2o_per_mton_n2o = c.ppt_per_mol / c.n2o_N_molar_mass * c.g_per_ton * c.ton_per_mton / c.ppt_per_ppb; 
c.time_const_for_n2o = 121;     % years
% Natural N2O emissions : AR6 WG1 Chapter 5. Table 5.3 Global N2O budget [ TgN / yr ]
% averaged over the 1980s, 1990s, 2000s as well as the recent decade starting in 2007.   
c.natural_n2o_emissions = 11.2;  % Mton / year

% initial_n2o_conc = 308.725;
c.initial_n2o_in_atm = c.N0 / c.ppb_n2o_per_mton_n2o; %  [ Mtons ]


% CH4 ............................................................
c.a3      =   -8.9603e-5;       % W / m^2 / ppb
c.b3      =   -1.2462e-4;       % W / m^2 / ppb
c.d3      =    0.045194;        % W / m^2 / ppb^(0.5)
c.M0      = 1714.0;             % initial CH4 concentration - based on 1990 levels from GISS, 2020. [ ppm ]  
c.ch4_ref =  731.41;            % reference CH4 concentration [ ppb ]
c.ch4_molar_mass = 16;          % g / mole
c.ppb_ch4_per_mton_ch4 = c.ppt_per_mol / c.ch4_molar_mass * c.g_per_ton * c.ton_per_mton / c.ppt_per_ppb; 

% Time Const for CH4 [ years ]
% IPCC AR6. WG1. Table 6.2
% Methane lifetime due to chemical losses, soil uptake and total atmospheric lifetime based on CMIP6 multi-model analysis,
% and bottom-up and top-down methane budget estimates in Table 5.2.
% https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-6/
c.time_const_for_ch4 = 10.3;

% Initial CH4 in atm [ Mtons ]
c.initial_ch4_in_atm = c.M0 / c.ppb_ch4_per_mton_ch4;

% ODS ............................................................
% Approximate radiative forcing from aggregate averaged total gases regulated by Montreal Protocol: CFCs, HCFCs etc.
% Parameters for ODS averaged from AR6 table 7.SM.7 https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf.
c.ods_molar_mass  = 120;        % g / mole

% PFC ............................................................
% IPCC AR6. WG1. Table 7.SM.7. https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf
% Preindustrial pfc conc : pfc concentration in atmosphere prior to 1750
c.preindustrial_pfc_conc    = 40 ;   % ppt
c.pfc_radiative_efficiency  = 0.099; % watt/(ppb*meter*meter)

% SF6 ............................................................
% IPCC AR6. WG1. Table 7.SM.7. https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf
c.sf6_molar_mass = 146 ;             % g / mole  
% Preindustrial sf6 conc : sf6 concentration prior to 1750, assumed zero
c.preindustrial_sf6_conc = 0 ;       % ppt 
c.sf6_radiative_efficiency = 0.567 ; % watt/(ppb*meter*meter) 

% HFC ............................................................
% http://www.qc.ec.gc.ca/dpe/publication/enjeux_ges/hfc134a_a.html
% IPCC AR6. WG1. Table 7.SM.7. https://www.ipcc.ch/report/ar6/wg1/downloads/report/IPCC_AR6_WGI_FGD_Chapter07_SM.pdf
c.hfc_molar_mass = [ 102 , 70 , 52 , 120 , 84 , 66 , 170 , 134 , 252 ] ; % [ g/mole ]
% Preindustrial HFC conc : Concentration of each HFC in atmosphere pre-1750, assumed zero
c.preindustrial_hfc_conc = 0; % [ ppt ]
% HFC radiative efficiency
c.hfc_radiative_efficiency = [ 0.167, 0.191, 0.111, 0.234, 0.168, 0.102, 0.237, 0.24, 0.357 ] ; % [ watt/(ppb*meter*meter) ]

% H2 ............................................................
%Sand, M., Skeie, R.B., Sandstad, M. et al. A multi-model assessment of the Global Warming Potential of hydrogen. Commun Earth Environ 4, 203 (2023). https://doi.org/10.1038/s43247-023-00857-8
c.h2_storage_leakage_rate = 2; %		[ percent/year ]
% Average utilization of storage capacity : On average, storage capacity is utilized approximately 50% over the range of variation from 0 to maximum
c.average_utilization_of_storage_capacity  = 0.5; 
% Energy intensity of hydrogen : https://h2tools.org/hyarc/calculator-tools/hydrogen-heating-values-mass-basis
c.energy_intensity_of_hydrogen = 120; % [ MJ/kg H2 ]

c.ch4_erf_per_h2_flux = 0.00046;          %  [ watts/(meter*meter)/(tg h2/year) ]
c.o3_erf_per_h2_flux = 0.00040;           %  [ watts/(meter*meter)/(tg h2/year) ]
c.strat_h2o_erf_per_h2_flux = 0.00019;    %  [ watts/(meter*meter)/(tg h2/year) ]
c.aerosol_rf_per_h2_flux = 0;             %  [ watts/(meter*meter)/(tg h2/year) ]

% Initial C concentration atmosphere and ocean [ ppm ]
c.C_init = [ c.C0 ; 1044.4 ; 3114.46 ; 3099.9 ; 13361.1 ; 18483 ]; % Carbon in each layer

% Preind Ocean C per meter : Corresponds with 767.8 GtC in a 75m layer.  [ GtonsC / meter ]
c.preind_ocean_c_per_meter = 10.2373;

% Preind C in Mixed Layer : Initial carbon concentration of mixed ocean layer.  [ GtonsC ]
c.preind_c_in_mixed_layer = c.preind_ocean_c_per_meter * c.layer_depth(1);


% Sensitivity of C Uptake to Temperature : (dimensionless)
% Allows users to control the strength of the feedback effect of temperature on uptake of C by land and oceans.
% 0 means no temperature-carbon uptake feedback and default of 1 yields the average value found in
% Friedlingstein et al., 2006.
% Climate-Carbon Cycle Feedback Analysis: Results from the C4MIP Model Intercomparison. Journal of Climate. p3337-3353.
c.sensitivity_of_c_uptake_to_temperature = 1;

% Sensitivity of pCO2 DIC to Temperature Mean [ percent/DegreesC ]
% Sensitivity of equilibrium concentration of dissolved inorganic carbon to temperature.
% Calibrated to be consistent with Friedlingstein et al., 2006.
% Climate-Carbon Cycle Feedback Analysis: Results from the C4MIP Model Intercomparison.
% Journal of Climate. p3337-3353.
% Default Sensitivity of C Uptake to Temperature of 1 corresponds to mean value from the 11 models tested.
c.sensitivity_of_pco2_dic_to_temperature_mean = 0.3;

% Sensitivity of pCO2 DIC to Temperature ( 1 / DegreesC )
% Sensitivity of pCO2 of dissolved inorganic carbon in ocean to temperature.
c.sensitivity_of_pco2_dic_to_temperature = c.sensitivity_of_c_uptake_to_temperature * c.sensitivity_of_pco2_dic_to_temperature_mean / c.c100_percent;


% Flux Ocean C to CH4 as CH4 [ Mtons/Year ]
% Weber, T., Wiseman, N.A. & Kock, A.
% Global ocean methane emissions dominated by shallow coastal waters. Nat Commun 10, 4584 (2019).
% https://doi.org/10.1038/s41467-019-12541-7. The Global Methane Budget 2000-2017.
% Saunous et al. 2020.
% https://essd.copernicus.org/articles/12/1561/2020/. Proxy for ocean and freshwater emissions.
c.flux_ocean_c_to_ch4_as_ch4 = 30;

% Flux Ocean C to CH4 : Flux of Ocean C to CH4 converted from CH4 per year [ GtonsC/Year ]
c.flux_ocean_c_to_ch4 = c.flux_ocean_c_to_ch4_as_ch4 / c.ch4_per_c / c.mton_c_per_gton_c;

% Buff C Coeff : Coefficient of CO2 concentration influence on buffer factor.  (dimensionless)
c.buff_c_coeff = 3.92;

% Ref Buffer Factor : Normal buffer factor.  (dimensionless)
c.ref_buffer_factor = 9.7;

c.buffer_factor_init = c.ref_buffer_factor; 

% =======================================================================
%{

% Reference RTE : Percent round trip efficiency for which the relationship between coverage required and renewables share is determined.
reference_rte = 100; 

Component weight[Elec Paths,Component]  =	1,1,0,0,0,0,0,0,0;

Capacity adjustment time	8		years

						1,1,0,0,0,0,0,0,0;
						1,1,0,0,0,0,0,0,0;
						1,1,0,0,0,0,0,0,0;
						0,1,1,0,0,0,0,0,0;
						0,1,0,1,0,0,0,0,0;
						0,0,0,0,1,0,0,0,0;
						0,0,0,0,0,1,0,0,0;
						0,0,0,0,0,0,1,0,0;
						0,0,0,0,0,0,0,1,0;
						0,0,0,0,0,0,0,0,1; (dimensionless)

Construction capability adj time	5	 years
adjusts the time to complete construction to reflect capacitated delays due to the limits of construction materials/resources.

Constr time sensitivity 	 1		dimensionless
Sensitivity of how much longer construction takes as afunction of loading, where a value of 1 indicates a linear relationship.

Future Plant construction time for electricity [ years ]
ECoal: 6 , EOil: 4 , EGas: 2 , EBio: 4 , nuclear: 6 , Renewable types: 2 , new: 5 , hydro: 4 
The number of years between starting to construct the new resource unit until the time that construction is completed; assumes permitting, planning and financing precede construction. Reference and data in the excel file En-ROADS Supporting Data A.

Historic plant construction time for electricity[Elec Paths]  = 6,4,2,4,6,4,2,2,2,2,5 year
Plant construction times from historic data. Reference and data in the excel file En-ROADS Supporting Data A.

Max exp term			10		dimensionless

Min construction capability	1		(EJ/year)/year
Min unit revenue  		1		$/GJ

Normal construction utilization 0.7		(dimensionless)



Retirement effect max  		8		dimensionless
Retirement effect at zero PCCR  1.1		dimensionless
Retirement effect exponent  	2 		dimensionless


SL adjustment time		2		years
Strat h2O erf per h2 flux	0.00019		watts/(meter*meter)/(tg h2/year)
Supply line recognition		0.5		unitless

TandD costs  			0.02	 $/kWh
Costs for transmission and distribution per unit generation, assumed constant. Estimated from https://www.eia.gov/todayinenergy/detail.php?id=32812: 2.2 - 3.2 cents/kWh in 2016 Also: https://energy.utexas.edu/sites/default/files/UTAustin_FCe_TDA_2016.pdf: average distribution capital costs approximately 0.6–0.8 ¢/kWh and average distribution system O&M costs approximately 0.4 ¢/kWh; average transmission costs approxiamtely 0.3-0.9 ¢/kWh.
equations

%}
