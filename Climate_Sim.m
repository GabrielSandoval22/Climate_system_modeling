% Climate_Sim.m
% simulate the EnRoads climate module

EnRoads_Constants; % in a structure c.xyz

% initial values for the state variables -------------

x0 = [ c.heat_initial; % heat in atmosphere and ocean  layers  (surface)  Watt.year/sq.m
       c.C0 / c.ppm_co2_per_GtonC ;      %  C  mass Gton
       c.N0 / c.ppb_n2o_per_mton_n2o  ;  % N2O mass Mton 
       c.M0 / c.ppb_ch4_per_mton_ch4 ];  % CH4 mass Mton

t = [ 1990 : 1.0 : 2100 ];          % time, years
N = length(t);                      % number of time points
u = zeros(1,N);                     % exogenous inputs

% solve the system of ode's -----------------------

[ t , x , dtxt , y ] = ode4u('Climate_Sys' , t , x0 , u , c );

% outputs ------------------------------------------

temperatures = y;

% plots --------------------------------------------

formatPlot(18,3,3) 
figure(1)
 clf
 semilogy(t,x)
 ylabel('states')
 legend('Q_a', 'Q_2', 'Q_3', 'Q_4', 'Q_5', 'C_a', 'N_2O', 'CH_4' )

figure(2)
 clf
 plot(t,y)
 ylabel('temperature outputs')
 legend('T_a', 'T_2', 'T_3', 'T_4', 'T_5' );

% ======================================================= 
