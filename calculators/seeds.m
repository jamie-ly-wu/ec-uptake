clear;
clc;
% This is a script for calculating the number of seeds and reagent volumes
% needed.

% Input date of experiment
date = '200621';

%% Calculations
% Input number of particles needed
np_num = 1e15;
% Input safety factor
sf = 1.5;
% Total number of particles to synthesize
np_to_synth = np_num*sf;

% Number of np from 100 mL of seeds synthesis
seeds_in_100ml = 1.45e14;
% Output total volume of seeds required (unit = mL)
vol_seeds_total = np_to_synth/seeds_in_100ml*100;

% Surface area of one 15 nm np
SA_per_np = 4*pi*(15/2)^2;
% Total surface area
SA_total = SA_per_np*np_to_synth;
% Input PEG length (mw)
peg_mw = 5000;
% Input PEG number per nm^2
peg_per_nm2 = 4;
% Number of PEG molecules needed
peg_num = SA_total*peg_per_nm2;
% Mass of PEG needed
peg_mass = peg_num/(6.022e23)*peg_mw;

%% Export to file
fileID = fopen('np_synthesis_calc.txt','w');

fprintf(fileID, 'Synthesis of seeds on %s \n\n', date);
fprintf(fileID, 'To make %1.3e seeds with a safety factor of %1.1f:\n', np_num, sf);
fprintf(fileID,'Total volume of seeds required is: %4.2f (mL)\n', vol_seeds_total);
fprintf(fileID,'Total mass of PEG-%5.0f required is: %2.4f (g)\n', peg_mw, peg_mass);

fclose(fileID);
