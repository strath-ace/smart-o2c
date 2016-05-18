function mjd2000 = date2mjd2000(date)

%DATE2MJD2000 Modified Julian day 2000 number from Gregorian date.
%
% MJD2000 = DATE2MJD2000(DATE) returns the modified Julian day 2000 number
% of the given date (Gregorian calendar) plus a fractional part depending
% on the time of day.
%
% Note: The function is valid for the whole range of dates since 12:00 noon
%       24 November -4713, Gregorian calendar. (This bound is set in order
%       to have symmetry with the inverse function JD2DATE)
% Note: The inputs must be feasible (i.e. the date must exist!). If an
%       unfeasible date is inputed, wrong results are given because no
%       check is done on that.
%
% INPUT
%	date        Date in the Gregorian calendar, as a 6-element vector
%               [year, month, day, hour, minute, second]. For dates before
%               1582, the resulting date components are valid only in the
%               Gregorian proleptic calendar. This is based on the
%               Gregorian calendar but extended to cover dates before its
%               introduction. date must be after 12:00 noon, 24 November
%               -4713.
%
% OUTPUT
%   mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% See also MJD20002DATE.
%
% FUNCTIONS CALLED
%    DATE2JD
%
% - Nicolas Croisard - 16/02/2008
% - Revised by Camilla Colombo - 03/03/2008
%
% ------------------------- - SpaceART Toolbox - --------------------------

mjd2000 = date2jd(date) - 2400000.5 - 51544.5;


return