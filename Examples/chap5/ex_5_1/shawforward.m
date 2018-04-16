% Calculate the Shaw forward problem for a model m
% for data points in s and equidistant model points in th    
% from Parameter Estimation and Inverse Problems, 3rd edition, 2013
% by R. Aster, B. Borchers, C. Thurber
%
    function g = shawforward(s,mbasis,th,dth)
    m=length(s);
    g=zeros(m,1);
    for i=1:m
        z=sin(s(i))+sin(th);
            q=sin(pi*z)./(pi*z);
            %patch here for the sin(x)/x limiting points
            q(z==0)=1;
            g(i) = sum(((cos(s(i))+cos(th)).^2).*(q.^2).*mbasis)*dth;   
    end
    
