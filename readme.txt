This is the code and data I used for: 

Palmer, G. 2017, An input-output based net-energy assessment of an electricity supply industry, Energy

The two Matlab scripts are I_O.mat and load_data.m.

There is data for several years. To select the year of interest, change the first line of 'load_data.m' to '2013-14' etc. Only run one year at a time.

There are several options in the main file 'I_O.mat'

Line 24 has the 'mode' default is '0' for all fuels, but can be changed. See lines 28 to 43.

Line 45 has 'CED', options are 0 or 1. Default is 0, which treats electricity as a fuel and used for EROI. Option 1 processess the electricity further into the primary feedstock fuels combusted for electricity generation, and gives the CED for fossil fuels. For calculation of CED for renewable fuels, see: 
Palmer and Floyd, 2017, An Exploration of Divergence in EPBT and EROI for Solar Photovoltaics, https://doi.org/10.1007/s41247-017-0033-0

Line 198 has the industry sector. The script can run for a single, multiple, or list of sectors. Use usual Matlab syntax.

Lines 514 to 617 write to file.

For queries, email me at graham.palmer@gmail.com