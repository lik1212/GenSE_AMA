%% Demostration file for general SE
%  As input you need z_all_data, z_all_flag, LineInfo and U_eva, check
%  the input data description for more details. 
%
% Author(s):    R. Brandalik
%
% Contact: brandalikrobert@gmail.com, brandalik@eit.uni-kl.de
%
% Special thanks go to the entire TUK ESEM team.
%
% Parts of the work were the result of the project CheapFlex, sponsored by
% the German Federal Ministry of Economic Affairs and Energy as part of the
% 6th Energy Research Programme of the German Federal Government. 

%% Clear start

path(pathdef); clear; close; clc

%% Path preperation

addpath([pwd,'\Subfunctions']);  % Add subfunction path

%% Load Demo Data

Input_Prep           = struct   ;
Input_Prep.Grid_Name = 'S1a_de' ;

% load([pwd,'\Demo_Data\Demo_Data_', Grid, '.mat']); 
load([pwd,'\Demo_Data\Demo_Data_', Input_Prep.Grid_Name, '_noisy.mat']); 

%% Inputs for State Estimation (can be extended with Inputs)

Inputs_SE.max_iter = 10         ; % Max num of iteration
Inputs_SE.z_conv   = 1 * 10^-0  ; % Abort criterion (convergence limit)
Inputs_SE.U_start  = 400/sqrt(3); % Voltage of iteration start (Flat-Start)

%% Reduce measurements

% time_steps = 1:100;
% z_all_data = z_all_data(:,time_steps);

%% In demo data change virtual measurement to real with some sigma

z_all_flag.Sigma     (z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type ~= 2) = 1;
z_all_flag.Accur_Type(z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type ~= 2) = 1;

%% Main estimation

tic
[x_hat, z_hat, z_hat_full, Out_Optional] = GenSE(z_all_data, z_all_flag, LineInfo, Inputs_SE);
toc