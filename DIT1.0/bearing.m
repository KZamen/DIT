function ang=bearing(start_lat,start_lon,end_lat,end_lon)
	%Calculates the bearing between two geographic coordinates in degrees
	%bearing=bearing(start_lat,start_lon,end_lat,end_lon)
	% check: http://gis.stackexchange.com/questions/29239/calculate-bearing-between-two-decimal-gps-coordinates
    % test code below on http://www.sunearthtools.com/tools/distance.php

	
	startlat=deg2radians(start_lat);
	startlon=deg2radians(start_lon);
	endlat=deg2radians(end_lat);
	endlon=deg2radians(end_lon);
	
	dlon=endlon-startlon;
	dphi=log(tan(endlat/2.0+pi/4.0)/tan(startlat/2.0+pi/4.0));
	
	if abs(dlon)>pi
		if (dlon>0.0)
			dlon= -(2.0*pi - pi);
		else
			dlon=(2.0*pi+dlon);
		end
	end
	bearing=mod(((atan2(dlon,dphi)+360)*180/pi),360.0);
	ang=deg2radians(bearing);
end