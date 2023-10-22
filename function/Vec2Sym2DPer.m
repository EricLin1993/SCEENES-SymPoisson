function [ X_Symetric ] = Vec2Sym2DPer( xc,c )
% global TC
% tic


uflag=triu(ones(c),0);
X_Symetric = zeros(c);

X_Symetric(logical(uflag)) = xc;

Temp = (X_Symetric+X_Symetric.');
X_Symetric = Temp-0.5*diag(diag(Temp));



%    it3=0;
%     for it2=1:c
%           for it1=1:it2
%                it3 = it3+1;
%                X_Symetric(it1,it2)= xc(it3);
%                if it1~=it2
%                    X_Symetric(it2,it1)= xc(it3);
%               end
%           end
%     end

% t =  toc;
% TC = TC+t;

end

