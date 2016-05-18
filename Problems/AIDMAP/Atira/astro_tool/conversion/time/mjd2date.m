function date = mjd2date(mjd)

%MJD2DATE Gregorian calendar date from modified Julian day number.
%
% DATE = MJD2DATE(MJD) returns the Gregorian calendar date (year, month,
% day, hour, minute, and second) corresponding to the modified Julian day
% number.
%
% Note: Since this function calls jd2date, MJD cannot be less than 
%       -2400000.5, that is 24 November -4713, 12:00 noon, Gregorian
%       proleptic calendar.
%
% INPUT
%   mjd         Date in modified Julian Day. The MJD count is from 00:00
%               midnight at the beginning of Wednesday November 17, 1858.
%               It must be a real greater or equal than -2400000.5.
%
% OUTPUT
%   date        Date in the Gregorian calendar, as a 6-element vector
%               [year, month, day, hour, minute, second]. For dates before
%               1582, the resulting date components are valid only in the
%               Gregorian proleptic calendar. This is based on the
%               Gregorian calendar but extended to cover dates before its
%               introduction.
%
% See also DATE2MJD.
%
% FUNCTIONS CALLED
%   MJD2JD
%   JD2DATE
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 21/02/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------


% Compute the year month and day
jd   = mjd2jd(mjd);
date = jd2date(jd);


return