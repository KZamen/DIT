function [mL1,p_value]=IRLS(G,d,tol)
    % Computing L1 regression IRLS (Iteratively reweighted Least Squares)
    % [model,P]=IRLS(G,d,tol)
    % model= Obtained model using L1 regression
    % P= misfit measure
    % G= Coeficient Matrix
    % d= Data Vector
    % tol= Tolerance
    % 
    % Kamyar Zamen May 8 2015
    
    eps=10E-5;
    
    m0=(G'*G)\(G'*d);
    res=d-(G*m0);
    R=zeros(size(res,1));
    
    for ii=1:size(res,1)
        if res(ii,1) < eps
            R(ii,ii) = 1/eps;
        else
            R(ii,ii)=res(ii,1);
        end
    end
    
    mL1_1=m0;
    mL1_2=(G'*R*G)\(G'*R*d);
    
    box=(norm((mL1_2-mL1_1),2))/(1+norm(mL1_2,2));
    
    while (box<tol)
        mL1_1=mL1_2;
        res=d-(G*mL1_1);
        R=zeros(size(res,1));
    
        for ii=1:size(res,1)
            if res(ii,1) < eps
                R(ii,ii) = 1/eps;
            else
                R(ii,ii)=res(ii,1);
            end
        end
        mL1_2=(G'*R*G)\(G'*R*d);
    end
    
    mL1=mL1_2;
    
    %%%% Calculating P-Value
    
    v=size(G,1)-size(G,2);
    sigma1=sqrt((1-(2/pi))*v);
    mue=norm((G*mL1-d),1);
    x=(mue-(sqrt(2/pi)*v))/sigma1;
    gama=(2-(pi/2))/((pi/2-1)^(3/2)*sqrt(v));
    S=normcdf(x,0,sigma1);
    Z=(((x^2)-1)/sqrt(2*pi))*exp(-(x^2)/2);
    p_value=1-S+(gama*Z/6);

end