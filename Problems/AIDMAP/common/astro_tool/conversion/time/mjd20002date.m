function date = mjd20002date(mjd2000)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
% MJD20002DATE Gregorian calendar date from modified Julian day 2000 number.
% 
% DATE = MJD20002DATE(MJD2000) returns the Gregorian calendar date
% (year, month, day, hour, minute, and second) corresponding to the
% modified Julian day 2000 number.
% 
% Note: Since this function calls jd2date, MJD2000 cannot be less than
%      -2451545, that is 24 November -4713, 12:00 noon, Gregorian
%      proleptic calendar.
% 
% INPUT
% 	mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%              since 01-01-2000, 12:00 noon. It must be a real greater or
%              equal than -2451545.
% 
% OUTPUT
%  date        Date in the Gregorian calendar, as a 6-element vector
%              [year, month, day, hour, minute, second]. For dates before
%              1582, the resulting date components are valid only in the
%              Gregorian proleptic calendar. This is based on the
%              Gregorian calendar but extended to cover dates before its
%              introduction.
% 
% See also DATE2MJD2000.
% 
% FUNCTIONS CALLED
%  MJD20002JD
%  JD2DATE
% 
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 21/02/2008
% 
% ------------------------- - SpaceART Toolbox - --------------------------

% Compute the year month and day
jd   = mjd20002jd(mjd2000);
date = jd2date(jd);


return
