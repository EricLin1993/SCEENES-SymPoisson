function [v,k] = Poisson1D(s,p,z)
%s=1.0 p=64 z=256
ld = z / p;
w = 2.0;


while 1
    i = 1;
    k = 1;
    
    while i < z
        v(k) = i;
        i = i+1;
        k = k+1;
        M_PI = 3.14;
        i = i+poisson((ld-1.0)*w*sin(M_PI*i/(z+1)));
    end   
    k = k-1;
    
    if(k > p)
        w = w*1.02;
    end
    if(k < p)
        w = w/1.02;
    end    
    
    if(k == p)
        break;
    end 
    clear v
end
% for k = 0:p-1
%     fprintf("%d\n",v(k)); 
% end
% v
end

%% subfunction
function [k] = poisson(lambd)
L = exp(-lambd); k = 0; p = 1;
while 1
    u = rand;
    p = p * u;
    k = k + 1;
    if(p < L)
        break;
    end
end
k = k - 1;
end

