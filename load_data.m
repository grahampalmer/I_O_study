dirName = '2013-14';
load(strcat(dirName,'/Z.mat'));  % our main intermediate table
load(strcat(dirName,'/y.mat'));  % final demand
load(strcat(dirName,'/ABS_4604.mat')); % energy consumption for BREE sectors
load(strcat(dirName,'/electricity_feedstock.mat')); % electricity supply
load(strcat(dirName,'/BREE_IOIG.mat'));  % mapping table between National and Energy Accounts
load(strcat(dirName,'/ANZIC_size.mat'));  % size of the National Accounts table
load(strcat(dirName,'/INVENTORY.mat')); % need to adjust final demand with inventory
load(strcat(dirName,'/ANZIC_names.mat'));    % ANZIC names
load(strcat(dirName,'/ANZIC_names_short.mat'));  % shortened names
load(strcat(dirName,'/BREE_types.mat')); % BREE type names
load('primary_energy_multiplier.mat');  % multipliers for each of 15 fuel types
load('fuel_types.mat'); % fuel type names
