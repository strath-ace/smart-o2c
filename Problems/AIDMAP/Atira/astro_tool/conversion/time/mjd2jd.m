function jd = mjd2jd(mjd)

%MJD2JD Modified Julian day number from Julian day number.
%
% JD = MJD2JD(MJD) returns the Julian day number corresponding to the given
% modified Julian day number.
%
% INPUT
%   mjd         Date in modified Julian Day. The MJD count is from 00:00
%               midnight at the beginning of Wednesday November 17, 1858.
%
% OUTPUT
%   jd          Date in Julian Day. The JD (Julian day) count is from 0 at
%               12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%               calendar. The corresponding date in Gregorian calendar is
%               12:00 noon, 24 November -4713.
%
% See also JD2MJD.
%
% FUNCTIONS CALLED
%   none
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 20/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


jd = mjd + 2400000.5;


return