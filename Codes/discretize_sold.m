function [ sold ] = discretize_sold( a )
%sold=[0.4,0.5,0.6,0.8,1,1.2,1.4,1.5,1.6,1.8,2];

if a>=0.4&&a<=0.5
    lower_diff=abs(a-0.4);
    upper_diff=abs(a-0.5);
    if lower_diff<=upper_diff
        sold=0.4;
    else
        sold=0.5;
    end
elseif a>0.5&&a<=0.6
    lower_diff=abs(a-0.5);
    upper_diff=abs(a-0.6);
    if lower_diff<=upper_diff
        sold=0.5;
    else
        sold=0.6;
    end
elseif a>0.6&&a<=0.8
    lower_diff=abs(a-0.6);
    upper_diff=abs(a-0.8);
    if lower_diff<=upper_diff
        sold=0.6;
    else
        sold=0.8;
    end
elseif a>0.8&&a<=1
    lower_diff=abs(a-0.8);
    upper_diff=abs(a-1);
    if lower_diff<=upper_diff
        sold=0.8;
    else
        sold=1;
    end
elseif a>1&&a<=1.2
    lower_diff=abs(a-1);
    upper_diff=abs(a-1.2);
    if lower_diff<=upper_diff
        sold=1;
    else
        sold=1.2;
    end
elseif a>1.2&&a<=1.4
    lower_diff=abs(a-1.2);
    upper_diff=abs(a-1.4);
    if lower_diff<=upper_diff
        sold=1.2;
    else
        sold=1.4;
    end
elseif a>1.4&&a<=1.5
    lower_diff=abs(a-1.4);
    upper_diff=abs(a-1.5);
    if lower_diff<=upper_diff
        sold=1.4;
    else
        sold=1.5;
    end
elseif a>1.5&&a<=1.6
    lower_diff=abs(a-1.5);
    upper_diff=abs(a-1.6);
    if lower_diff<=upper_diff
        sold=1.5;
    else
        sold=1.6;
    end
elseif a>1.6&&a<=1.8
    lower_diff=abs(a-1.6);
    upper_diff=abs(a-1.8);
    if lower_diff<=upper_diff
        sold=1.6;
    else
        sold=1.8;
    end
elseif a>1.8&&a<=2
    lower_diff=abs(a-1.8);
    upper_diff=abs(a-2);
    if lower_diff<=upper_diff
        sold=1.8;
    else
        sold=2;
    end
    
end


end

