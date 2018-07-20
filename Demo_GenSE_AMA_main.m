%% Demostration file for general SE with the Augmented Matrix approach
%  As input you need z_all_data, z_all_flag, LineInfo and U_eva, check
%  the input data description for more details. 
%
% Author(s):    R. Brandalik
%               D. Henschel
%               J. Tu
%               P. Gassler
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

load([pwd,'\Demo_Data\Demo_Data_v2.mat']); 

%% Inputs for State Estimation (can be extended with Inputs)

Inputs_SE.max_iter = 10         ; % Max num of iteration
Inputs_SE.z_conv   = 1 * 10^-1  ; % Abort criterion (convergence limit)
Inputs_SE.U_start  = 400/sqrt(3); % Voltage of iteration start (Flat-Start)

%% Main estimation

tic
[x_hat, z_hat, z_hat_full, Out_Optional] = GenSE_AMA(z_all_data, z_all_flag, LineInfo, Inputs_SE);
toc