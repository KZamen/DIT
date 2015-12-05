function distance=Haversine(lat1,lon1,lat2,lon2)
    %Calculates the distance between two Coordinates (in degrees)
    % []=Haversine(Latitude1, Longitude1,Latitude2, Longitude2)
	%check: http://www.movable-type.co.uk/scripts/latlong.html
    % test code below on http://www.sunearthtools.com/tools/distance.php

    
    a=6384;
    b=6353;
    avglat=(lat1+lat2)/2;
    R=sqrt((((a^4).*(cosd(avglat)).^2)+((b^4).*(sind(avglat)).^2))/((a.*cosd(avglat)).^2+(b.*sind(avglat)).^2));  % Earth radius at the desired latitude
    
    latd=lat2-lat1;
    lond=lon2-lon1;
    d=(sind(latd./2)).^2 + cosd(lat1) .* cosd(lat2) .* (sind(lond./2)).^2;
    distance=2.*R.*asind(d);
    
end