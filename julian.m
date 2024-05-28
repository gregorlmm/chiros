% FUNCTION TO CALCULATE THE JULIAN DAY NUMBER AND THE JULIAN DATE
% Valid (at least) for all non-negative Julian Day Numbers
% * i.e. since noon of 1 January 4713 BC in the proleptic Julian calendar;
% * i.e. since noon of 24 November 4714 BC in the proleptic Gregorian calendar;
% * i.e. since '24-Nov--4713 12:00:00' in the proleptic ISO 8601 calendar, Universal Time (UT1)
% INPUT: days
%   * as a MATLAB serial date number, i.e. days since '0-Jan-0000 00:00:00' or as a DateString; 
%   * Date should be in the proleptic ISO 8601 calendar;
%   * Time should be in Universal Time (UT1)
% OUTPUT: [JDN,JD]
%   * Julian Day Number
%   * Julian Date

function [JDN,JD] = julian(days)

	datevector = datevec(days);
	a = floor((14-datevector(2))/12);
	y = datevector(1)+4800-a;
	m = datevector(2)+12*a-3;
	JDN = datevector(3) + floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100)+floor(y/400)-32045;
	JD = JDN + (datevector(4)-12)/24 + datevector(5)/1440 + datevector(6)/86400;

end

% Algorithm taken from http://www.cs.utsa.edu/~cs1063/projects/Spring2011/Project1/jdn-explanation.html
% The first day for MATLAB is that which starts on 31 December 2 BC (proleptic Gregorian calendar)
% This date can be expressed as '0-Jan-0000 00:00:00', proleptic ISO 8601 calendar, Universal Time (UT1)...
% which is equivalent to the date '31-Dec--0001 00:00:00' (proleptic ISO 8601 calendar, UT1)