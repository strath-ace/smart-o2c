function mjd = mjd20002mjd(mjd2000)

%MJD20002JD Modified Julian day number from modified Julian day 2000
%number.
%
% MJD = MJD20002MJD(MJD2000) returns the modified Julian day number
% corresponding to the given modified Julian day 2000 number.
%
% INPUT
%	mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% OUTPUT
%   mjd         Date in modified Julian Day. The MJD count is from 00:00
%               midnight at the beginning of Wednesday November 17, 1858.
%
% See also MJD2MJD2000.
%
% FUNCTIONS CALLED
%   none
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 20/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


mjd = mjd2000 + 51544.5;


return