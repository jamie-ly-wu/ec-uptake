clear;
clc;
% This script calculates the number of cells to seed.

% Input date of experiment
date = '200621';

%% Calculations

% Doubling times for different cell types:
%doubling_times = py.dict(pyargs('b16f1_tyr-1-',18,'htb14',40.5));
% Input cell type
celltype = 'b16f1_tyr-1-';
time_doubling = 18;%doubling_times{celltype};
% Input target cell number at harvest
target_num = 6e6;
% Input safety factor
sf = 1;
% Input time (hr) between seeding and harvest
time_lapsed = 18;
% Output number of cells to seed
seed_num = target_num*sf/(2^(time_lapsed/time_doubling));


%% Export to file
fileID = fopen('cell_seeding_calc.txt','w');

fprintf(fileID, 'Cell seeding of %s on %s \n\n', celltype, date);
fprintf(fileID, 'To get %1.2e cells with a safety factor of %1.1f after %3.1f hours:\n', target_num, sf, time_lapsed);
fprintf(fileID,'Seed %1.2e cells\n', seed_num);

fclose(fileID);