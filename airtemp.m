function [T_air] = airtemp(T_amp,omega,delta,T_0,t_instant)

	T_air = T_0 + T_amp*sin(day(datetime(datestr(t_instant)),'dayofyear')/omega - delta);

end