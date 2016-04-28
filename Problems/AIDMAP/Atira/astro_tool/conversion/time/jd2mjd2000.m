function mjd2000 = jd2mjd2000(jd)

%JD2MJD2000 Modified Julian day 2000 number from Julian day number.
%
% [MJD2000] = JD2MJD2000(JD) returns the modified Julian day 2000 number
% corresponding to the given Julian day number.
%
% INPUT
%   jd          Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% OUTPUT
%   mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% See also MJD20002JD.
%
% FUNCTIONS CALLED
%   none
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Camilla Colombo - 29/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


mjd     =  jd - 2400000.5;
mjd2000 = mjd - 51544.5;


return