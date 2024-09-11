function [ mask ] = SymMask( m,upf,c)
mask =zeros(size(c));
idx=1;
for i = 1 : c
    for j = 1 : i
        if(m(idx) == 1)
            if(upf(idx) == 1)
                mask(j,i) = 1;
            else
                mask(i,j) = 1;
            end
        end
        idx=idx+1;
    end
end
end