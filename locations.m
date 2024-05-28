function [locparams] = locations(loc)
	locparams = table;

	switch(loc)
		case 1	% STOCKHOLM
			locparams.locname = "Stockholm";
			
            % AIR TEMPERATURE AND RADIATION MODEL
            % Geographical coordinates and time zone
			locparams.lat = 59.33; % degrees (verified)
			locparams.long = 18.07; % degrees (verified)
			locparams.t_zone = 1; % UTC+1 (verified)
			% Climate
			%locparams.T_min = -1; % Celsius
            %locparams.T_min = -1.8;
            locparams.T_min = -2.5; % Celsius (verified)
            %locparams.T_max = 18; % Celsius
            locparams.T_max = 18.2; % Celsius (verified)
			locparams.T_0 = (locparams.T_max + locparams.T_min)/2; % Celsius (verified)
			locparams.T_amp = (locparams.T_max - locparams.T_min)/2; % Celsius (verified)
			locparams.omega = 365/(2*pi); % 58.091554228541803 days (reported as 58.092, verified)
			locparams.doy_peak = 200; % DOY of peak summer day (verified)
			locparams.delta = 2*pi/365*locparams.doy_peak - pi/2; % 1.872044937413096 (reported as 1.872, verified)
			
            T_min = locparams.T_min;
            T_max = locparams.T_max;
            doy_peak = locparams.doy_peak;

            % DIAPAUSE MODEL
            % Base temperature for local chironomids (°C)
			%locparams.T_diap = 4;
            locparams.T_diap = 7.2; % verified
			% Light requirement to avoid diapause (minutes)
			%locparams.daylength_diap = 9*60;
            locparams.daylength_diap = 8*60; % 6 hours = 480 min, verified

		case 2	% BERLIN
			locparams.locname = "Berlin";

            % AIR TEMPERATURE AND RADIATION MODEL
			% Geographical coordinates and time zone
			locparams.lat = 52.52; % verified
			locparams.long = 13.405; % verified
			locparams.t_zone = 1; % UTC+1 (verified)
			% Climate
			locparams.T_min = 0.5; % Celsius (verified)
			locparams.T_max = 19.5; % Celsius (verified)
			locparams.T_0 = (locparams.T_max + locparams.T_min)/2; % Celsius (verified)
			locparams.T_amp = (locparams.T_max - locparams.T_min)/2; % Celsius (verified)
			locparams.omega = 365/(2*pi); % 58.091554228541803 days (reported as 58.092, verified)
			locparams.doy_peak = 200; % DOY of peak summer day (verified)
			locparams.delta = 2*pi/365*locparams.doy_peak - pi/2; % 1.872044937413096 (reported as 1.872, verified)

            T_min = locparams.T_min;
            T_max = locparams.T_max;
            doy_peak = locparams.doy_peak;

            % DIAPAUSE MODEL
			% Base temperature for local chironomids (°C)
			locparams.T_diap = 7.5; % verified
			% Light requirement to avoid diapause (minutes)
			locparams.daylength_diap = 9*60; % 9 hours = 540 min, verified

		case 3	% TRENTO
			locparams.locname = "Trento";

            % AIR TEMPERATURE AND RADIATION MODEL
			% Geographical coordinates and time zone
			locparams.lat = 46.0748; % reported as 46.075 (verified)
			locparams.long = 11.1217; % reported as 11.122 (verified)
			locparams.t_zone = 1; % UTC+1 (verified)
			% Climate
			locparams.T_min = 2; % Celsius (verified)
            locparams.T_max = 23.2; % Celsius (verified)
			locparams.T_0 = (locparams.T_min + locparams.T_max)/2; 
			locparams.T_amp = (locparams.T_max - locparams.T_min)/2;
			locparams.omega = 365/(2*pi);
			locparams.doy_peak = 200; % DOY of peak summer day (verified)
			locparams.delta = 2*pi/365*locparams.doy_peak - pi/2;
			
            T_min = locparams.T_min;
            T_max = locparams.T_max;
            doy_peak = locparams.doy_peak;

            % DIAPAUSE MODEL
            % Base temperature for local chironomids (°C)
			locparams.T_diap = 8;
			% Light requirement to avoid diapause (minutes)
			locparams.daylength_diap = 10*60;  % 10 hours = 600 min, verified		
	end
end