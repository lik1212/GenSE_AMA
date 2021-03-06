function Settings = defaultSettings(Inputs)
% defaultSettings checks the content of Inputs and based on it define the
% Settings. If some optional field do not occur in "Inputs" set the default
% values for it.
%
% Author(s): R. Brandalik

%% Transfer main info from "Inputs

Settings = Inputs;

%% Default setting if the inputs do not occur in the variable "Inputs"

default_Options = {...
    'max_iter' , 10          ;...
    'z_conv'   , 5 * 10^-6   ;...
    'U_start'  , 400/sqrt(3) ;...
    };

for k_Opt = 1 : size(default_Options,1)
    if ~isfield(Inputs, default_Options{k_Opt,1})
        Settings.(default_Options{k_Opt,1}) = default_Options{k_Opt,2};
    end
end

