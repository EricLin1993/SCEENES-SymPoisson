function [ res ] = Sym2DToVecPer( xc )
% global TC
% tic

uflag=triu(ones(size(xc)),0);
flag=uflag-diag(0.5*ones(size(xc,1),1));
X = (xc+xc.').*flag;

res = X(logical(uflag));

%    it3=0;
%     for it2=1:size(xc,2)
%           for it1=1:it2
%                it3 = it3+1;              
%                if it1~=it2
%                    res(it3) = (xc(it1,it2)+xc(it2,it1));
%                else    
%                    res(it3) = xc(it1,it2);
%                end
%           end
%     end
    

    
% t =  toc;
% TC = TC+t;
end

