This is the code and data I used for: 

Palmer, G. 2017, An input-output based net-energy assessment of an electricity supply industry, Energy

The two Matlab scripts are I_O.mat and load_data.m.

There is data for several years. To select the year of interest, change the first line of 'load_data.m' to '2013-14' etc. Only run one year at a time.

There are several options in the main file 'I_O.mat'

Line 24 has the 'mode'm default is '0' for all fuels, but can be changed. See lines 28 to 43.

Line 45 has 'CED', options are 0 or 1. Default is 0, which treats electricity as a fuel. Option 1 processess the electricity further into the primary feedstock fuels combusted for electricity generation.

Line 198 has the industry sector. The script can run for a single, multiple, or list of sectors. Use usual Matlab syntax.

Lines 514 to 617 write to file.

For queries, email me at graham.palmer@gmail.com