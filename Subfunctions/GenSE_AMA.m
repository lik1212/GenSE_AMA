function [x_hat, z_hat, z_hat_full, Out_Optional] = GenSE_AMA(z_all_data, z_all_flag, LineInfo, Inputs_SE)
%GenSE_AMA General State Estimation with augmented matrix approach 
%   Input
%       z_all_data - The measurement data
%       z_all_flag - Basic information about the measurements (what is
%                    measured, where is it measured ...). 
%                    See Description_Of_Inputs
%       LineInfo   - Basic information about the lines (branches) of the
%                    grid
%       Inputs_SE  - .max_iter  - Max num of iteration
%                    .z_conv    - Abort criterion (convergence limit)
%                    .U_start   - Voltage of iteration start (Flat-Start)
%
%   Output
%       x_hat        - Estimated state vector
%       z_hat        - Estimated measurement vector
%       z_hat_full   - Estimation of all important measurements:
%                      (U, phi, P & Q) in this order for all node in ID
%                      order
%       Out_Optional - Optional output.
%
% Author(s): R. Brandalik
%
% Contact: brandalikrobert@gmail.com, brandalik@eit.uni-kl.de
%
% Special thanks go to the entire TUK ESEM team.
%
% Parts of the work were the result of the project CheapFlex, sponsored by
% the German Federal Ministry of Economic Affairs and Energy as part of the
% 6th Energy Research Programme of the German Federal Government. 

%% Inputs to Settings

if nargin < 4; Inputs_SE = struct; end
Settings = defaultSettings(Inputs_SE);
max_iter = Settings.max_iter;
z_conv   = Settings.z_conv  ;
U_start  = Settings.U_start ;

%% Sort measurement values in this order: real, pseudo, virtual

[~, z_order     ] = sort(z_all_flag.Accur_Type);
[~, z_order_back] = sort(z_order);               % Flag for sorting back
z_all_flag_order  = z_all_flag(z_order,:);
z_all_data_order  = z_all_data(z_order,:);

%% Prepare Y_L1L2L3 & H_index

[Y_012, Y_012_Node_ID] = LineInfo2Y_012(LineInfo);           % Admittance matrix in symmetrical components
Y_L1L2L3               = Y_012_to_Y_L1L2L3(Y_012);           % Admittance matrix in phase sequences
num_Nodes              = size(Y_L1L2L3,1)/3;                 % Get number of grid nodes
H_index = table;
H_index.Node1_ID  = repmat(Y_012_Node_ID,4,1);            	 % H_index contains the assumed position of 
H_index.Phase     = repmat ([1; 2; 3]   , 4 * num_Nodes, 1); % all meas. functions for U, phi, P and Q
H_index.Meas_Type = repelem([1; 2; 3; 4], 3 * num_Nodes, 1);

HC_flag = NaN(size(z_all_flag_order, 1), 1);                 % Initial vector to sort H_index in the same order as z
for k_z = 1 : numel(HC_flag)                                 % Create vector to sort H_index in the same order as z
    HC_flag(k_z) = find(...
        H_index.Node1_ID  == z_all_flag_order.Node1_ID (k_z) & ...
        H_index.Phase     == z_all_flag_order.Phase    (k_z) & ...
        H_index.Meas_Type == z_all_flag_order.Meas_Type(k_z));
end

HC_switch = find(z_all_flag_order.Accur_Type == 3, 1);          % HC_switch defines the position were HC_flag seperates virtual from real-pseudo values
H_flag    = HC_flag(1         : HC_switch - 1, :);              % Position of real-pseudo values in H in relation to z
C_flag    = HC_flag(HC_switch :           end, :);              % Position of virtual     values in H in relation to z
R         = diag(z_all_flag_order.Sigma(1 : HC_switch - 1).^2); % Covariance matrix 
alpha     = min(diag(R));

%% Initial output (results)

num_inst = size(z_all_data, 2); % Number of instances
% Flat start for state vector for all instances
x_hat = repmat([...
    repmat(U_start,3 * num_Nodes,1); ... % Voltage
    repmat([...                          % Angle
    z_all_data(z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type == 2 & z_all_flag.Phase == 1) ; ...
    z_all_data(z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type == 2 & z_all_flag.Phase == 2) ; ...
    z_all_data(z_all_flag.Accur_Type == 3 & z_all_flag.Meas_Type == 2 & z_all_flag.Phase == 3)], num_Nodes, 1)],...
    1, num_inst);

z_hat_full = NaN(size(H_index,1), num_inst);
z_hat      = NaN(numel(HC_flag) , num_inst);
if nargout > 3
    flag_conv  = false(num_inst, 1); % Flag for information if converged
    num_iter   = NaN  (num_inst, 1); % Number of iteration
end

%% SE solver

% Over all instances
for k_inst = 1 : num_inst
    x_k_hat = x_hat(:, k_inst); % initial state vector
    for k_iter = 1 : max_iter % iteration
        z_SE    = get_z_SE(Y_L1L2L3, Y_012_Node_ID, x_k_hat); % Get the z vector off all measurements of SE, not z_hat!
        z_k_hat = z_SE([H_flag; C_flag]);                     % Estimated measurement vector z_hat
        delta_z = z_all_data_order(:, k_inst) - z_k_hat;
        if all(abs(delta_z) < z_conv) % Convergence limit reached
            if nargout > 3            % For optional output
                flag_conv(k_inst) = true  ;
                num_iter (k_inst) = k_iter; 
            end
            break;
        end
%         sum(delta_z.^2)   % Temp
        H_SE = get_H_SE(Y_L1L2L3, Y_012_Node_ID, x_k_hat); % Get the H matrix of SE
        H_AM = H_SE(H_flag,:);                             % Matrix of real-pseudo values, AM stands for augmented matrix
        C_AM = H_SE(C_flag,:);                             % Matrix of virtual     values
        % Build Hachtel matrix of the augmented matrix approach
        Hachtel     = sparse([ ...
            1/alpha.*R                       , H_AM                             , zeros(size(H_AM,1),size(C_AM,1));...
            H_AM'                            , zeros(size(H_AM,2),size(H_AM,2)) , C_AM'                           ;...
            zeros(size(C_AM,1),size(H_AM,1)) , C_AM                             , zeros(size(C_AM,1),size(C_AM,1)) ...
            ]);
        rhs_AM                  = [      ...      % Create right hand side of the LSE for the augmented matrix approach
            delta_z(1 : HC_switch - 1)  ;...
            zeros(size(C_AM,2),1)       ;...
            delta_z(HC_switch :   end)  ;...
            ];
        solution_vector = Hachtel \ rhs_AM;
        delta_x(:,1)    = solution_vector(HC_switch : end - numel(C_flag)); % State vector correction
        x_k_hat = x_k_hat + delta_x; % New state vector
    end
    if nargout > 3 && flag_conv(k_inst) == false % For optional output
        num_iter (k_inst) = max_iter;
    end
    x_hat     (:, k_inst) = x_k_hat;                % Estimated state vector
    z_hat     (:, k_inst) = z_k_hat(z_order_back);  % Estimated measurement vector
    z_hat_full(:, k_inst) = z_SE;                   % Estimation of all important measurements
end

%% Optional output

if nargout > 3
    Out_Optional.flag_conv = flag_conv;
    Out_Optional.num_iter  = num_iter ;
    Out_Optional.Y_L1L2L3  = Y_L1L2L3;
end

