function mjd = jd2mjd(jd)

%JD2MJD Modified Julian day number from Julian day number.
%
% MJD = JD2MJD(JD) returns the modified Julian day number corresponding to
% the given Julian day number.
%
% INPUT
%   jd          Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% OUTPUT
%   mjd         Date in modified Julian Day. The MJD count is from 00:00
%               midnight at the beginning of Wednesday November 17, 1858.
%
% See also MJD2JD.
%
% FUNCTIONS CALLED
%   none
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Camilla Colombo - 29/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


mjd = jd - 2400000.5;


return