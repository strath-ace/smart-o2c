function mjd2000 = mjd2mjd2000(mjd)

%MJD2MJD2000 Modified Julian day 2000 number from Modified Julian day
%number.
%
% [MJD2000] = MJD2MJD2000(MJD) returns the modified Julian day 2000 number
% corresponding to the given modified Julian day number.
%
% INPUT
%   mjd         Date in modified Julian Day. The MJD count is from midnight
%               at the beginning of Wednesday November 17, 1858.
%
% OUTPUT
%   mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% See also MJD20002MJD.
%
% FUNCTIONS CALLED
%   none
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Camilla Colombo - 29/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


mjd2000 = mjd - 51544.5;


return