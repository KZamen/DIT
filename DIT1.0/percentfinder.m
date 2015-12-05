function y=percentfinder(percent)
%% This function finds the coresponding Z in normal distribution 
    check=true;
    for ii=(0:0.01:3.9)
        nor1=normcdf(-1*ii,0,1);
        nor2=normcdf(ii,0,1);
        ptile=nor2-nor1;
        ptile=roundtoper(ptile,3);
        if (ptile==percent)
            y=ii;
            check=false;
        end
    end
    if (check)
        disp('No Exact Match for this percentile')
    end
end