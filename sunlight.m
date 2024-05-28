function [sunlight_duration] = sunlight(lat,long,datenum_UT1)

	[Julian_Day_Number,Julian_Date] = julian(datenum_UT1);
	Julian_Century = (Julian_Date-2451545)/36525;
	Geom_Mean_Long_Sun_deg = mod(280.46646 + Julian_Century*(36000.76983 + Julian_Century*0.0003032),360);
	Geom_Mean_Anom_Sun_deg = 357.52911 + Julian_Century*(35999.05029 - 0.0001537*Julian_Century);
	Sun_Eq_Ctr = sin(deg2rad(Geom_Mean_Anom_Sun_deg))*(1.914602 - Julian_Century*(0.004817+0.000014*Julian_Century)) + ...
                             sin(deg2rad(2*Geom_Mean_Anom_Sun_deg))*(0.019993 - 0.000101*Julian_Century) + ...
                             sin(deg2rad(3*Geom_Mean_Anom_Sun_deg))*0.000289;
	Sun_True_Long_deg = Geom_Mean_Long_Sun_deg + Sun_Eq_Ctr;
	Sun_App_Long_deg = Sun_True_Long_deg - 0.00569 - 0.00478*sin(deg2rad(125.04-1934.136*Julian_Century));
	Mean_Obliq_Ecliptic_deg = 23 + (26 + ((21.448 - Julian_Century*(46.815 + Julian_Century*(0.00059 - Julian_Century*0.001813))))/60)/60;
	Obliq_Corr_deg = Mean_Obliq_Ecliptic_deg + 0.00256*cos(deg2rad(125.04-1934.136*Julian_Century));
	Sun_Declin_deg = rad2deg(asin(sin(deg2rad(Obliq_Corr_deg))*sin(deg2rad(Sun_App_Long_deg))));
	HA_sunrise_deg = rad2deg(acos(cos(deg2rad(90.833))/(cos(deg2rad(lat))*cos(deg2rad(Sun_Declin_deg))) - tan(deg2rad(lat))*tan(deg2rad(Sun_Declin_deg))));
	sunlight_duration = 8*HA_sunrise_deg;

end