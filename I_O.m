%% Clear, load data, and allocate memory to matrices

clear all;

format long g;

doWaitBar = 0;  % 1 = TRUE 0 = FALSE

disp(['starting ...']);

load_data;  % load the data matrices and vectors

direct_total = 0;
n_ = ANZIC_size;    % number of sectors in ABS 5209 National Accounts
f_ = 15;    % number of BREE/ABS4604 fuel types
b_ = 24;    % number of BREE/ABS4604 sectors
t_ = 6;     % tiers + 1
mode = 0;   % default to zero, define some modes, 0: EROI, 1, 2, 3 .. respective fuels, 30: liquid fuels, 40: coal

pathway_cumulative_matrix = zeros(n_, f_);   % this holds each fuel type for the last sector run
                                                    % need to declare
                                                    % before mode loop so that it doesn't get
                                                    % cleared each time
pathway_tier_matrix = zeros(t_, f_);
                                                    
for mode = [10]     % run through every fuel type. normally 1:15

switch mode
        case 0  % mode 0 is default, all fuels
            disp('all fuels');
        case num2cell(1:f_) % all individual fuels
            disp(fuel_types(mode));
        case 30 % liquid fuels
            disp('... liquid fuels => 7,8,9,10,11');
        case 40 % black and brown coal, coke, coal by-products, briquettes
            disp('... black and brown coal, coke, coal by-products, briquettes => 1,2,3,4,5');
        otherwise
            fprintf('_mode_ not valid. check line 9');
            if doWaitBar
                close (h);
            end
            return;
end

CED = 0;    % 0=EROI, 1=cumulative energy demand, include primary energy for electricity
threshold = 0.0001;  % this is PJ, 0.001 = 1 TJ
output_rows = 1000000;
if doWaitBar
    h=waitbar(0, 'Working ...');
end
x = sum(Z, 2) + y; % total output = intermediate + final demand
                                % need to work out what to do with negative
                                % inventory change
v = x' - sum(Z,1);
% value added vector
% value added equals total output - sum of intermediate use column for each
% sector
                                
A = Z./repmat(x',n_,1);
% Aij = Zij / yj , total requirements matrix A
% ./ apply division element-wise
% repmat : repeat copies of array, transpose X
L = inv(eye(n_) - A);  % Leontief inverse
fuel_allocation_vector = zeros(n_,1);    % declare a vector for fuel intensity
fuel_allocation_matrix = zeros(n_,f_);   % this saves the individual fuels for each sector

summary = zeros(n_,8);             % summary table for each industry for each tier
clear F;
e_out = zeros(1, n_);

F = zeros(n_, n_);          % this is the output matrix for each fuel type
F_total = zeros(n_, 144);   % this is the cumulative total for all fuel types
output_matrix_size = 150000; % set at 150000
pathway_ANZIC = cell (output_matrix_size, 1);
tier_0_total = 0;
tier_1_total = 0;
tier_2_total = 0;
tier_3_total = 0;
tier_4_total = 0;
tier_5_total = 0;
count = int64(1);   % use a 64 bit integer to economise memory

numbered_vector = cell(output_matrix_size, 1);
% fill a vector with consecutive numbers
for i = 1:output_matrix_size
    numbered_vector{i} = sprintf('%d', i);
end

e = zeros(1, n_);  % fuel per unit of output

%% Allocate BREE energy data to ABS industry sectors

for fuel = drange(1:f_) % we have 15 fuel types in ABS 4609/BREE
    
    for BREE_sector = drange(1:b_)
        % BREE broad industry sectors 'agriculture' 'mining' 'food beverages alcohol'
        % etc ... #22 is net use by housholds
        total = 0;
        % sum the columns that are part of the BREE sector and divide into the
        % quantity of fuel for that sector -> PJ per $
        for i = drange(1:n_)   % ABS 5209 industry sectors
            if BREE_IOIG(i) == BREE_sector
                % total = total + sum(Z(:,i));    % Z is our 114x114 transactions matrix
                total = total + x(i);  % changed Jan 2017
                % sum the total supplies for each of the sectors within the
                % Energy Account group
                
            end
            
        end
                       
        for i = drange(1:n_)
            if BREE_IOIG(i) == BREE_sector
                % each of the 114 ABS industry sectors gets the same fuel
                % intensity as a proportion of that industry relative to
                % all within the BREE industry group
                % change Dec 2016, was sum(Z(:,i))
                % intensity_vector(i) = (sum(Z(:,i)) / total) * ABS_4604_table(BREE_sector, fuel);
                fuel_allocation_vector(i) = x(i) / total * ABS_4604_table(BREE_sector, fuel);
                fuel_allocation_matrix(i, fuel) = fuel_allocation_vector(i);
                % each sector (114) rows, and each fuel type (15)
            end
        end
    end
end

%% Re-distribute primary energy for electricity generation if CED = 1

if CED == 1
    % if CED mode, needing to extract electricity feedstock and add to
    % vectors, and set electricity to zero
    for i = 1:n_
        if i ~= 65 % if not electricity generation sector
            (sum(fuel_allocation_matrix(i,15))/sum(fuel_allocation_matrix(:,15)) * electricity_feedstock) + fuel_allocation_matrix(i,:);
            % set the proportion of electricity for this sector relative to
            % total, multiplied by the electricity feedstock vector, then
            % add back to the intensity vector
            fuel_allocation_matrix(i,15) = 0; % now set electricity to zero
        else
            fuel_allocation_matrix(65,:) = fuel_allocation_matrix(65,:) + electricity_feedstock;
            % but for electricity generation, just set to feedstock
            fuel_allocation_matrix(65,15) = 0;
            % and make electricity zero
        end
    end
end

%% Add together fuel types depending on mode, apply primary energy multiplier and allocate to vector 'e'

for fuel = drange(1:f_) % we have 15 fuel types in ABS 4609/BREE
    
    switch mode
        case 0  % mode 0 is default, all fuels
            e_out = fuel_allocation_matrix(:,fuel)'*inv(diag(x)); % e_out is a 1x114 vector fuel per unit of total output
        case num2cell(1:f_) % all individual fuels
            if fuel == mode; e_out = fuel_allocation_matrix(:,fuel)'*inv(diag(x));
            else e_out = 0;
            end
        case 30 % liquid fuels
            if fuel == 7 || fuel == 8 || fuel == 9 || fuel == 10 || fuel == 11; e_out = fuel_allocation_matrix(:,fuel)'*inv(diag(x));
            else e_out = 0;
            end
        case 40 % black and brown coal, coke, coal by-products, briquettes
            if fuel == 1 || fuel == 2 || fuel == 3 || fuel == 4 || fuel == 5; e_out = fuel_allocation_matrix(:,fuel)'*inv(diag(x));
            else e_out = 0;
            end
        otherwise
            fprintf('_mode_ not valid. check line 9');
            if doWaitBar
                close (h);
            end
            return;
    end
    
    
    if mode == 0    % keep adding fuels as loop runs through all fuels
        e = primary_energy_multiplier(fuel) * e_out + e;
    else
        e = e_out + e;  % don't use multiplier for specific fuels
    end
    
    % e is a 1x114 vector for total fuel for each industry
    % accumulate each e, multiply by fuel multiplier (coal = 1,
    % elec = 3, etc, we have 15 fuel types
    
    % F = diag(e_out)*L*diag(y);   % output matrix - use of fuel
    
end

if doWaitBar
    close(h);
end
if doWaitBar
    h = waitbar(0, 'Working ...');
end
%% Create energy pathway data and table for each ABS sector

for sector = [66]    % pick some sectors or do the lot, electtricity is [72] for 109 sectors, 65/66 for 111 & 114 sectors
    
    tic;    % start a timer
    
    if doWaitBar
        waitbar(0);
    end
    
    if sector > n_
        fprintf('Sector not valid. Stopping');
        close (h);
        return;
    end
    
    switch mode % save the energy for the given mode to use to calculate GVA intensity
        case 0  % mode 0 is default -> all fuels
            mode_energy = sum(fuel_allocation_matrix(sector,:)); % all fuels
        case num2cell(1:f_)
            mode_energy = sum(fuel_allocation_matrix(sector,mode));    % just the specific fuel
        case 30 % liquid fuels
            mode_energy = sum(sum(fuel_allocation_matrix(sector,7:11))); % fuels 7 thru 11
        case 40 % black and brown coal
            mode_energy = sum(sum(fuel_allocation_matrix(sector,1:2))); % fuels 1 and 2
    end
    
    % allocate some memory for large matrices
    pathway_energy = zeros (output_matrix_size, 1); % this is the stored value
    pathway_cumulative_vector = zeros (n_, 1); % cumulate for each sector for tier 1
    pathway_text = cell  (output_matrix_size, 1);   % this is the pathway stored as text
    pathway_ANZIC = cell (output_matrix_size, 1);
    pathway_cumulative_text = cell (output_matrix_size, 1);  % tier 1 has 114 sectors
    count = 1;
    % fprintf('\n');
    disp([sector, ANZIC_names(sector)]);    % output to command line
    
    % tier 0 is direct, we're running this one sector at a time
    % ///// this was e x y //////
    temp = e(sector) * y(sector);
    % direct energy = fuel per unit of output(e) x final demand (y)
    % the rest of the energy attributed to that sector becomes embedded in the intermediate use matrix (Z) 
    
    tier_0 = temp;
    disp(['direct ',num2str(fuel_allocation_matrix(sector,mode))]);
    direct_total = direct_total + fuel_allocation_matrix(sector,mode);
    % y is final demand
    
    tier_1 = 0;
    tier_2 = 0;
    tier_3 = 0;
    tier_4 = 0;
    tier_5 = 0;
    % tier 1 is taking the direct energy inputs of every sector that
    % supplies the target sector
    
    tier_0_total = tier_0_total + temp;
    pathway_cumulative_vector (sector) = tier_0_total;  % start with direct (tier_0)
    pathway_text {count}  = 'Direct';   % store the pathway as text
    pathway_energy (count)  = temp;     % store the result
    count = count + 1;
    
    
    if doWaitBar
        close(h);
        h = waitbar(5 / 145, 'Extracting pathways ...');
    end    
    %% This is our pathway coding for 5 tiers
    % Start at the top tier and truncate below if less than the threshold
    for j = drange(1:n_)
        if doWaitBar
            waitbar(0.05 + j / 145);
        end
        
        if j == sector
            if j == 65 || j == 66
                tier_0 = tier_0/(y(j)/x(j));
                pathway_cumulative_vector (j) = tier_0;  % start with direct (tier_0)
                disp(['sector 65/66, tier_0 : ',num2str(tier_0)]);
                % for electricity, all of the direct energy is attributed
                % since tier_0 is already calculated, use the (Z/total)
            end
        end
        
        temp = e(j) * A(j, sector) * y(sector);
        % e : 1x114, A : 114x114, y : 114x1
        % so there are 114 paths for each sector at tier 1
        % e is fuel per unit of output
        
        if j == sector && (j == 65 || j == 66)
            temp = 0;
        end
        % since we are adding all of sectors 65 or 66, we need to
        % subtract these from tiers to avoid double counting
        
        
        if temp > threshold     % truncate below threshold
            tier_1 = tier_1 + temp;
            pathway_text {count}  = num2str(j); % store the pathway as text
            pathway_energy (count)  = temp;       % store the result
            pathway_cumulative_vector (j) = pathway_cumulative_vector (j) + temp;
            count = count + 1;
            for i = drange(1, n_)  % this is the inputs (rows) from each sector
                temp = e(i) * A(i, j) * A(j, sector) * y(sector);
                % i is the input sector
                
                if i == sector && (i == 65 || i == 66)
                    temp = 0;
                end
                % since we are adding all of sectors 65 or 66, we need to
                % subtract these from tiers to avoid double counting
                
                if temp > threshold
                    tier_2 = tier_2 + temp;
                    pathway_text {count} = [num2str(j), '-', num2str(i)];  % store the pathway as text
                    pathway_energy (count) = temp;   % store the result
                    pathway_cumulative_vector (j) = pathway_cumulative_vector (j) + temp;
                    count = count + 1;
                    for k = drange(1, n_)
                        temp = e(k) * A(k, i) * A(i, j) * A(j, sector) * y(sector);
                        
                        if k == sector && (k == 65 || k == 66)
                            temp = 0;
                        end
                        % since we are adding all of sectors 65 or 66, we need to
                        % subtract these from tiers to avoid double counting
                        
                        if temp > threshold
                            tier_3 = tier_3 + temp;
                            pathway_text {count} = [num2str(j), '-', num2str(i), '-', num2str(k)];  % store the pathway as text
                            pathway_energy (count) = temp;   % store the result
                            pathway_cumulative_vector (j) = pathway_cumulative_vector (j) + temp;
                            count = count + 1;
                            for l = drange(1, n_)
                                temp = e(l) * A(l, k) * A(k, i) * A(i, j) * A(j, sector) * y(sector);
                                
                                if l == sector && (l == 65 || l == 66)
                                    temp = 0;
                                end
                                % since we are adding all of sectors 65 or 66, we need to
                                % subtract these from tiers to avoid double counting
                                
                                if temp > threshold
                                    tier_4 = tier_4 + temp;
                                    pathway_text {count} = [num2str(j), '-', num2str(i), '-', num2str(k), '-', num2str(l)];  % store the pathway as text
                                    pathway_energy (count) = temp;   % store the result
                                    pathway_cumulative_vector (j) = pathway_cumulative_vector (j) + temp;
                                    count = count + 1;
                                    for m = drange(1, n_)
                                        temp = e(l) * A(m, l) * A(l, k) * A(k, i) * A(i, j) * A(j, sector) * y(sector);
                                        
                                        if m == sector && (m == 65 || m == 66)
                                            temp = 0;
                                        end
                                        % since we are adding all of sectors 65 or 66, we need to
                                        % subtract these from tiers to avoid double counting
                                        
                                        if temp > threshold
                                            tier_5 = tier_5 + temp;
                                            pathway_text {count} = [num2str(j), '-', num2str(i), '-', num2str(k), '-', num2str(l), '-', num2str(m)];  % store the pathway as text
                                            pathway_energy (count) = temp;   % store the result
                                            pathway_cumulative_vector (j) = pathway_cumulative_vector (j) + temp;
                                            count = count + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Aggregate tiers for each fuel type
    
    pathway_tier_matrix(1, mode) = tier_0;
    pathway_tier_matrix(2, mode) = tier_1;
    pathway_tier_matrix(3, mode) = tier_2;
    pathway_tier_matrix(4, mode) = tier_3;
    pathway_tier_matrix(5, mode) = tier_4;
    pathway_tier_matrix(6, mode) = tier_5;
     
    
    %% Fill tier 1 pathway totals for each fuel type, exclude for mode zero (all)
    
    if mode > 0
        pathway_cumulative_matrix(:, mode) = pathway_cumulative_vector;  % copy the entire column across
    end
    
    
    %% Output progress to command screen
    
    % output to command screen
    disp(['Tier 0 ', num2str(tier_0)]);
    tier_1_total = tier_1_total + tier_1;
    disp(['Tier 1 ', num2str(tier_1)]);
    tier_2_total = tier_2_total + tier_2;
    disp(['Tier 2 ', num2str(tier_2)]);
    tier_3_total = tier_3_total + tier_3;
    disp(['Tier 3 ', num2str(tier_3)]);
    tier_4_total = tier_4_total + tier_4;
    disp(['Tier 4 ', num2str(tier_4)]);
    tier_5_total = tier_5_total + tier_5;
    disp(['Tier 5 ', num2str(tier_5)]);
    disp(['Total ', num2str(tier_0+tier_1+tier_2+tier_3+tier_4+tier_5)]);
    summary(sector, 1) = sector;
    summary(sector, 2) = tier_0+tier_1+tier_2+tier_3+tier_4+tier_5;
    summary(sector, 3) = tier_0;
    summary(sector, 4) = tier_1;
    summary(sector, 5) = tier_2;
    summary(sector, 6) = tier_3;
    summary(sector, 7) = tier_4;
    summary(sector, 8) = tier_5;
    
    
    %% Calculate scaled tiers and total
    
    % these convert energy (PJ) to intensity (GJ/$m)
    % /// these were X, should they be y ???
    tier_0_scaled = 1000 * tier_0 / y(sector);  % divide by total $, multiply by 1000 to give GJ/$000
    tier_1_scaled = 1000 * tier_1 / y(sector);  % start with PJ and $m
    tier_2_scaled = 1000 * tier_2 / y(sector);
    tier_3_scaled = 1000 * tier_3 / y(sector);
    tier_4_scaled = 1000 * tier_4 / y(sector);
    tier_5_scaled = 1000 * tier_5 / y(sector);
    total = tier_0 + tier_1 + tier_2 + tier_3 + tier_4 + tier_5;
    total_scaled = 1000 * total / y(sector);
    
    %% Sort pathways
    
    % now sort the pathways largest to smallest with an index to line up
    % the text with the result
    
   
    if doWaitBar
        close(h);
        h = waitbar(120/145, 'Sorting ...');
    end
    [pathway_energy_sorted, I_] = sort (pathway_energy, 'descend'); % sort largest to smallest with an index
    pathway_text_sorted = pathway_text(I_); % order text vector by index
    
    %% Extract sectors from sorted list to print sector names in output table
    
    
    if doWaitBar
        close(h);
        h = waitbar(125 / 145, 'Extracting sector names ...');
    end
    
    % need to extract the sector names from the results for printing in the
    % final output table
    for i = 1:size(pathway_text_sorted);
        if not(sum(pathway_text_sorted{i}) > 0)  % break at end of list
            break;
        end
        
        % only need to extract names for displayed pathways
        if i > output_rows
            break;
        end
        
        % pathway_cumulative_tier_1_text{1} = [num2str(total)];
        % CHECK THIS - CHANGE BECAUSE FIRST ITEM HAS PATHWAY ENERGY =
        % TOTAL
        % pathway_cumulative_tier_1_text{1} = ANZIC_names_short{1};
        
        % otherwise continue
        tempString = strsplit(pathway_text_sorted{i},'-');   % split text string of path into sectors
        [nr,nc] = size(tempString); % nr = number of rows, nc = number of columns
        if nc == 1    % only one item, this is an entry for tier 1
            if not(isnan(str2double(tempString(1))))    % check whether str2double returns a valid argument
                pathway_ANZIC{i} = [ANZIC_names_short{str2double(tempString(1))}];    % return the ANZIC name
                if i == 1.23456 %i == 1 CHECK THIS - JUST PUT THIS IN TO TRY
                    pathway_cumulative_text{i} = [num2str(total)]; % this is the first row - 'direct'
                    % pathway_cumulative_tier_1_text{i} = ANZIC_names_short{1};
                else
                    % pathway_cumulative_tier_1_text{i} = ANZIC_names_short{1};
                    pathway_cumulative_text{i} = num2str(round(...
                        pathway_cumulative_vector(str2num(pathway_text_sorted{i})), 6));
                    % pathway_cumulative_tier_1_text{i} = strtrim(cellstr(num2str(round(...
                    %    pathway_cumulative_tier_1(str2num(pathway_text_sorted{i})), 6))));
                end
            end
        elseif nc == 2  % two items, this is an entry for tier 2
            pathway_ANZIC{i} = [ANZIC_names_short{str2double(tempString(1))}, '/', ...
                ANZIC_names_short{str2double(tempString(2))}];    % concatenate the two
            pathway_cumulative_text{i} = ['-'];
        else  % three or more
            pathway_ANZIC{i} = [ANZIC_names_short{str2double(tempString(1))}, '/', ...
                ANZIC_names_short{str2double(tempString(2))}, '/', ...
                ANZIC_names_short{str2double(tempString(3))}] ;    % concatenate three
            
            pathway_cumulative_text{i} = ['-'];
        end
    end
    
    %% Convert sorted lists to cell arrays for concatenation and writing to file
    
    
    if doWaitBar
        close(h);
        h = waitbar(130 / 145, 'Converting to cell arrays ...');
    end
    
    pathway_energy_sorted_c  = strtrim(cellstr(num2str(round(pathway_energy_sorted, 6))));
    % convert numbers to cell array so these can be concatenated
    % this is PJ
    pathway_intensity_sorted_c = strtrim(cellstr(num2str(round((1000 / y(sector)) * pathway_energy_sorted, 6))));
    % convert numbers to cell array
    % this is GJ/$000
    
    % now concatenate the call arrays together to give us some rows to print
    %for j = drange(1:n_)
    %    if
    output_table = [numbered_vector pathway_text_sorted pathway_energy_sorted_c pathway_cumulative_text pathway_ANZIC];
    %end
    
    %% Write to file
    
    
    if doWaitBar
        close(h);
        h = waitbar(135 / 145, 'Writing to file ...');
    end
    
    % write to file
%    if CED == 1
%       fileID = fopen(['data/', num2str(sector), '_', num2str(mode), '_CED.dat'], 'w') ; % use /data/ folder, append CED to file
%        disp(['writing ','data/', num2str(sector), '_', num2str(mode), '_CED.dat']);
%    else
%        fileID = fopen(['data/', num2str(sector), '_', num2str(mode), '_EROI.dat'], 'w') ; % use /data/ folder, append EROI to file
%        disp(['writing ','data/', num2str(sector), '_', num2str(mode), '_EROI.dat']);
%    end
    
    if CED == 1
        fileID = fopen([strcat(dirName,'/data/'), num2str(sector), '_', num2str(mode), '_CED.dat'], 'w') ; % use /data/ folder, append CED to file
        disp(['writing ',strcat(dirName,'/data/'), num2str(sector), '_', num2str(mode), '_CED.dat']);
    else
        fileID = fopen([strcat(dirName,'/data/'), num2str(sector), '_', num2str(mode), '_EROI.dat'], 'w') ; % use /data/ folder, append EROI to file
        disp(['writing ',strcat(dirName,'/data/'), num2str(sector), '_', num2str(mode), '_EROI.dat']);
    end    
    
    formatSpec = '%s \n'; fprintf(fileID,formatSpec,['[',num2str(sector),']',ANZIC_names{sector}]);
    formatSpec = '%s \n'; fprintf(fileID,formatSpec,strcat('[5209] Australian National Accounts: Input-Output Tables: ', dirName));
%   Primary energy multipliers labels for the file print
%     formatSpec = '%s \n'; fprintf(fileID,formatSpec, '-----------------------------------------------------------------');
%     formatSpec = '%s \n'; fprintf(fileID,formatSpec,['Primary energy multipliers :']);
%     
%     for fc = drange(1:f_)   % fc - fuel counter, print fuel types and multipliers
%         formatSpec = '%20s %8.2f \n'; fprintf(fileID,formatSpec, fuel_types{fc}, primary_energy_multiplier(fc));
%     end
    
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, '-----------------------------------------------------------------');
    
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, ['Fuel mode: ', mode_str(mode)]);
    
    if CED == 1
        CED_str = 'Cumulative energy demand';
    else
        CED_str = 'EROI - exclude electricity primary feedstock';
    end
    
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, ['CED or EROI mode?: ', CED_str]);
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, ['Threshold: ', num2str(threshold, '%10.8f')]);
    
    for row = 1:output_matrix_size % count the total pathways so we can print
        if not(sum(output_table{row, 2}) > 0); break; end
    end
    
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, ['Found: ', comma(row), ' pathways']);
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, ['Elapsed time: ', num2str(toc), ' secs']);
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, '-----------------------------------------------------------------');
    formatSpec = '%s \n'; fprintf(fileID,formatSpec,'Input-output table use table rows :');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'Total supply             $', comma(x(sector)), 'm');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'Total industry uses      $', comma(sum(Z(sector,:))), 'm');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'Final demand             $', comma(y(sector)), 'm');
    formatSpec = '%s \n';        fprintf(fileID,formatSpec,'Input-output table use table columns :');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'Total intermediate uses  $', comma(sum(Z(:,sector))), 'm');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'GVA (total - interm.)    $', comma(x(sector) - sum(Z(:,sector))), 'm');
    formatSpec = '%s %8s %s \n'; fprintf(fileID,formatSpec,'GVA Intensity             ',...
        comma((1000000 * mode_energy / (x(sector) - sum(Z(:,sector))))), 'GJ/$m GVA');
    formatSpec = '%s \n'; fprintf(fileID, formatSpec, ' ');
    formatSpec = '%s \n'; fprintf(fileID,formatSpec,'Calculated totals (not Leontief) :');
    formatSpec = '%-6s %12s  %12s %12s \n'; fprintf(fileID, formatSpec, 'Tier', 'PJ', 'GJ/$1000', 'Proportion');
    formatSpec = '%s \n'; fprintf(fileID, formatSpec, ' ');
    
    formatSpec = '%-12s %8.4f %10.4f \n'; fprintf(fileID,formatSpec,'Total', total, total_scaled);
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 0', tier_0, tier_0_scaled, 100 * (tier_0 / total), '%');
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 1', tier_1, tier_1_scaled, 100 * (tier_1 / total), '%');
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 2', tier_2, tier_2_scaled, 100 * (tier_2 / total), '%');
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 3', tier_3, tier_3_scaled, 100 * (tier_3 / total), '%');
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 4', tier_4, tier_4_scaled, 100 * (tier_4 / total), '%');
    formatSpec = '%-12s %8.4f %10.4f %10.2f %s \n'; fprintf(fileID,formatSpec,'Tier 5', tier_5, tier_5_scaled, 100 * (tier_5 / total), '%');
    formatSpec = '%s \n'; fprintf(fileID,formatSpec, '-----------------------------------------------------------------');
    
    if doWaitBar
        waitbar(140 / 145);
    end
    
    % this is our main pathway table
    formatSpec = '%-6s %-15s %-10s %-10s %s \n'; fprintf(fileID, formatSpec, 'Item', 'Pathway', 'Tier', 'Pathway', 'Description first three');
    formatSpec = '%29s %10s \n'; fprintf(fileID, formatSpec, 'Energy', 'Energy');
    % formatSpec = '%30s \n'; fprintf(fileID, formatSpec,'(GJ/$1000)');
    formatSpec = '%29s %10s \n'; fprintf(fileID, formatSpec,'(PJ)', '(PJ)');
    formatSpec = '%-6s %-15s %-10s %-10s %s \n';
    
    [nrows,ncols] = size(output_table);
    if nrows > output_rows
        nrows = output_rows;
    end
    
    for row = 1:nrows  % print our data points
        if row > 100; break; end % only print first 100
        % if not(sum(output_table{row, 2}) > 0); break; end % only write up to threshold
        fprintf(fileID,formatSpec,output_table{row,:});
    end
    
    % close the file and end
    fclose(fileID) ;
    
end

%% Finish off with command-line summary

% these are the Leontief totals
tier_total = diag(e)*L*diag(y);
tier_0 = diag(e)*eye(n_)*diag(y);
tier_1 = e * diag(A * y);
tier_2 = e * diag(A * A * y);
tier_3 = e * diag(A * A * A * y);
tier_4 = diag(e)*A*A*A*A*diag(y);
tier_5 = diag(e)*A*A*A*A*A*diag(y);

% print to the command screen
disp([' ']);
fprintf ('Leontief totals\n');
disp(['Total ', num2str(sum(sum(tier_total))), '   ',...
    num2str(tier_0_total+tier_1_total+tier_2_total+tier_3_total+tier_4_total+tier_5_total)]);

disp(['Tier 0 ', num2str(sum(sum(tier_0))), '  ',num2str(tier_0_total)]);
disp(['Tier 1 ', num2str(sum(sum(tier_1))), '  ',num2str(tier_1_total)]);
disp(['Tier 2 ', num2str(sum(sum(tier_2))), '  ',num2str(tier_2_total)]);
disp(['Tier 3 ', num2str(sum(sum(tier_3))), '  ',num2str(tier_3_total)]);
disp(['Tier 4 ', num2str(sum(sum(tier_4))), '  ',num2str(tier_4_total)]);
disp(['Tier 5 ', num2str(sum(sum(tier_5))), '  ',num2str(tier_5_total)]);
disp(['Direct total ',num2str(direct_total)]);

if doWaitBar
    close (h);
end

end

dlmwrite(strcat(dirName,'/results_',num2str(sector),'.csv'),pathway_cumulative_matrix);  % write the sector-fuel matrix to a csv file
dlmwrite(strcat(dirName,'/tiers_',num2str(sector),'.csv'),pathway_tier_matrix);  % write the sector-fuel matrix to a csv file


