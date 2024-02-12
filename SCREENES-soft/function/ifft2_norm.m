function K = ifft2_norm(M)
% Two dimensional orthogonal IFFT
% y=myifft2(x)

[x,y] = size(M);
K = ifft2(ifftshift(ifftshift(M,1),2)) * sqrt(x) * sqrt(y);

end