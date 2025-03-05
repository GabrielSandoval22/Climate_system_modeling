% Carbon_Sim.m
% simulate the EnRoads Carbon Cycle module

EnRoads_Constants; % in a structure c.xyz

% initial values for the state variables -------------

x0 = [ c.C_init / c.ppm_co2_per_GtonC ;  %  C  mass Gton C
       c.N0 / c.ppb_n2o_per_mton_n2o  ;  % N2O mass Mton 
       c.M0 / c.ppb_ch4_per_mton_ch4  ]; % CH4 mass Mton

t = [ 1990 : 1.0 : 2100 ];          % time, years
N = length(t);                      % number of time points
u = temperatures;                   % exogenous inputs from Climate_Sim

% solve the system of ode's -----------------------

[ t , x , dtxt , y ] = ode4u('Carbon_Sys' , t , x0 , u , c );

% outputs ------------------------------------------

% nothing yet ...

% plots --------------------------------------------

formatPlot(18,3,3);

figure(3)
 clf
 plot(t,x)
 ylabel('states')
 legend('C_a', 'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'N_2O', 'CH_4')

% ======================================================= 
