%% Demostration file for general SE with compare
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

path(pathdef); clear; close all; clc

%% Path preperation

addpath([pwd,'\Subfunctions'        ]);  % Add subfunction path
addpath([pwd,'\Comparison_Functions']);  % Add comparison subfunction path

%% Input Prepareration

Input_Prep                  = struct   ;
Input_Prep.Grid_Name        = 'S1a_de' ;
Input_Prep.ResDate          = ''       ;
Input_Prep.LF_Res_Path      = [pwd,'\Comparison_Data\'];
Input_Prep.SinInfo_Path     = [pwd,'\Comparison_Data\'];
Input_Prep.SE_Inputs_Path   = [pwd,'\Demo_Data\'];
Input_Prep.NodeRes_Name     = [Input_Prep.Grid_Name, '_NodeRes_raw',   Input_Prep.ResDate, '.mat'];
Input_Prep.BranchRes_Name   = [Input_Prep.Grid_Name, '_BranchRes_raw', Input_Prep.ResDate, '.mat'];
Input_Prep.with_TR    = true                                                                     ;

if Input_Prep.with_TR
    Input_Prep.NodeRes_Name       = [Input_Prep.      NodeRes_Name(1 : end - 4) ,'_wo_TR.mat'];
    Input_Prep.BranchRes_Name     = [Input_Prep.    BranchRes_Name(1 : end - 4) ,'_wo_TR.mat'];
end

%% Load Demo Data

load([Input_Prep.SE_Inputs_Path, 'Demo_Data_', Input_Prep.Grid_Name, Input_Prep.ResDate,          '.mat']); 
load([Input_Prep.SE_Inputs_Path, 'Demo_Data_', Input_Prep.Grid_Name, Input_Prep.ResDate, '_noisy.','mat']); 

% Just for comparison
load([Input_Prep.SinInfo_Path, 'SinInfo_', Input_Prep.Grid_Name, '.mat']); 
NodeRes_all_exakt   = load([Input_Prep.LF_Res_Path, Input_Prep.NodeRes_Name  ]);
BranchRes_all_exakt = load([Input_Prep.LF_Res_Path, Input_Prep.BranchRes_Name]);

%% Inputs for State Estimation (can be extended with Inputs)

Inputs_SE.max_iter = 10         ; % Max num of iteration
Inputs_SE.z_conv   = 1 * 10^-0  ; % Abort criterion (convergence limit)
Inputs_SE.U_start  = 400/sqrt(3); % Voltage of iteration start (Flat-Start)

%% If reduced measurements time steps wanted

% time_steps = 550:569;
% z_all_data = z_all_data(:,time_steps);

%% In demo data change virtual measurement to real with some sigma

z_all_flag.Sigma     (z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type ~= 2) = 1;
z_all_flag.Accur_Type(z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type ~= 2) = 1;

%% Main estimation

tic
[x_hat, z_hat, z_hat_full, Out_Optional] = GenSE(z_all_data, z_all_flag, LineInfo, Inputs_SE);       % without AMA
% [x_hat, z_hat, z_hat_full, Out_Optional] = GenSE_AMA(z_all_data, z_all_flag, LineInfo, Inputs_SE); % with    AMA
toc

%% For comparison

NodeRes_all_estim   = z_full2NodeRes_all(z_hat_full, SinInfo);
BranchRes_all_estim = NodeRes2BranchRes(NodeRes_all_estim, SinInfo, Out_Optional.Y_L1L2L3);

NodeRes_all_exakt   = NodeRes_all_exakt  .NodeRes_all;
BranchRes_all_exakt = BranchRes_all_exakt.BranchRes_all;

if ~exist('time_steps', 'var') % If reduced time steps
    time_steps = unique(NodeRes_all_estim.ResTime)';
end

plot_comparison(NodeRes_all_estim, BranchRes_all_estim, NodeRes_all_exakt, BranchRes_all_exakt, time_steps);
