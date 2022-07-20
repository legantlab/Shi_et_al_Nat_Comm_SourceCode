load('Det_PSF_OTF_510_NA1p0_RichardsWolf.mat');
load('plane wave sets for GS Hex and SQ.mat');

NA_max = 0.21;               % maximum NA of the annulus
NA_min = 0.0;               % minimum NA of the annulus
NA_ideal = (NA_min+NA_max)/2;
NA_det = 1;                  % NA of the detection objective

xy_pol = [1 0];              % Polarizaiton of simulated beam

fill_factor = 1;             % Fill factor for the annulus, set to 1 the beam will occupy the entire annulus
crop_factor = 0.00;          % Crop factor on the SLM, regions with intensity less than the croping factor 
                             % relative to the peak intensity wil be
                             % cropped out. Not used for Gaussian beam
                             % simulations.

ny_step = 12;                % number of steps to simulate along the propagation direction
y_stepsize = 4;              % size of each step, unit in lambda

lattice_descrip = 'Gaussian';
PW =PW_Gaussian;

detPSF = xz_PSF_RW_510nm_NA1p0;
detOTF = xz_OTF_RW_510nm_NA1p0;

gamma = 0.5;                 % gamma factor used for plotting

folder = '.\simulation_output\Test';

%generate XZ profile along different positions in propagation
Calc_and_Plot_3D_LLS_PSFs_and_OTFs_Overall(lattice_descrip, xy_pol, PW, NA_max, NA_ideal, NA_min, NA_det, fill_factor, crop_factor, ny_step, y_stepsize, detPSF, detOTF, gamma,folder,1);
close all;


%generate YZ profile
ny_step = 501;                  % number of steps to simulate along the propagation direction
y_stepsize = 0.1;               % size of each step, unit in lambda

Calc_and_Plot_3D_LLS_PSFs_and_OTFs_Overall_YZprofile(lattice_descrip, xy_pol, PW, NA_max, NA_ideal, NA_min, NA_det, fill_factor, crop_factor, ny_step, y_stepsize, detPSF, detOTF, gamma,folder,1);
close all;