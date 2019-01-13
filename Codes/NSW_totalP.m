function [P_02_div_P_01_NSW]=NSW_totalP (M1)
if (M1<=1)
    P_02_div_P_01_NSW=1;
else 
    gamma=1.4;
    M2=sqrt((1+(((gamma-1)./2).*M1.^2))./((gamma.*M1.^2)-((gamma-1)./2)));
    a=(1+((gamma-1)./2).*M1.^2).^(gamma./(gamma-1));
    b=(1+((gamma-1)./2).*M2.^2).^(gamma./(gamma-1));
    P_02_div_P_01_NSW=(b./a).*(1+((2.*gamma)./(gamma+1)).*((M1.^2)-1));
end
end

