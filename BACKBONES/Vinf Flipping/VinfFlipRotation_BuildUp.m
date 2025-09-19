function [Vinf_Flip_Seqs, Vinf_Flip_FlybysData] = VinfFlipRotation_BuildUp(VinfFlip_TransfersDatabase, PseudoReso_TransfersDatabase, pars)

% This function computes a petal rotation sequence given a set of flyby
% options and the desired change in argument of periapsis.

% Author: J.C. Garcia Mateas
% Last revision: 29/05/2024

%% INPUTS %%
%
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %%
%
% - Vinf_Flip_Seqs: 

%% CHANGELOG %%


%% BUILD THE PETAL ROTATION SEQUENCES %%
% Initialize Flag to indicate when desired petal rotation is achieved
Flag_PetalRot_BuildUP = 0;

% Initialize Flag regarding first segment of the Vinf Flip Sequence
Flag_Flyby_1 = 1;

% Initialize Flyby Counter
Flyby_Count = 0;


while Flag_PetalRot_BuildUP == 0

    if Flag_Flyby_1 == 1 % If the first flyby must be computed

        % Retrieve Resonant & Pseudo-Resonant Transfers to be considered
        Reso_Options    = VinfFlip_TransfersDatabase;
        PseudoReso_Opts =  PseudoReso_TransfersDatabase;

        % Initialize variables
        Buffer_FlybyData = [];

        tic
        % Explore possible flybys & construct petal rotation sequences
        [Buffer_Transfers, Buffer_Transfers_IDX, Buffer_FlybyData] = Transfer_Combinations(Reso_Options, PseudoReso_Opts, Buffer_FlybyData, pars);

        % Rank solutions according to change in true longitude
        [Buffer_Transfers, Buffer_FlybyData] = Ranking_Solutions(Buffer_Transfers, Buffer_FlybyData, pars);

        toc

        % Change Flag of First join after accommplished
        Flag_Flyby_1 = 0;

        Flyby_Count = Flyby_Count + 1;

         % Check if maximum number of flybys has been reached
         if Flyby_Count == pars.PetalRot.Num_Flybys 
             Flag_PetalRot_BuildUP = 1;
         end

    else

        % Convert the Transfers Buffer into 3D Matrix Form 
        Buff_Transfers_Matrix = cat(3, Buffer_Transfers{:});

        tic

        if mod(Flyby_Count, 2) ~= 1   
            % A pseudo-resonance must be considered
            [Buffer_Transfers, Buffer_Transfers_IDX, Buffer_FlybyData] = Transfer_Combinations(Buff_Transfers_Matrix, PseudoReso_Opts, Buffer_FlybyData, pars);

            % Rank solutions according to change in true longitude
            [Buffer_Transfers, Buffer_FlybyData] = Ranking_Solutions(Buffer_Transfers, Buffer_FlybyData, pars);

        else 
            % A Vinf Flip resonance must be considered
            [Buffer_Transfers, Buffer_Transfers_IDX, Buffer_FlybyData] = Transfer_Combinations(Buff_Transfers_Matrix, Reso_Options, Buffer_FlybyData, pars);

            % Rank solutions according to change in true longitude
            [Buffer_Transfers, Buffer_FlybyData] = Ranking_Solutions(Buffer_Transfers, Buffer_FlybyData, pars);
        end

        toc

         Flyby_Count = Flyby_Count + 1;

         % Check if maximum number of flybys has been reached
         if Flyby_Count == pars.PetalRot.Num_Flybys 
             Flag_PetalRot_BuildUP = 1;
         end

    end

end


%% PREPARE THE FUNCTION OUTPUT %%
% Define function outputs
Vinf_Flip_Seqs       = Buffer_Transfers;
Vinf_Flip_FlybysData = Buffer_FlybyData;

end