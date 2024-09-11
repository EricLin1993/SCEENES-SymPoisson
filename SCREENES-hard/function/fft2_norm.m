function K = fft2_norm(M)
% Two dimensional orthogonal FFT
% y=myfft2(x)

[x,y] = size(M);
K = fftshift(fftshift(fft2(M),1),2) / sqrt(x) / sqrt(y);

end