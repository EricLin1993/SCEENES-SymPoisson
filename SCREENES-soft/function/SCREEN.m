function Rec_X = SCREEN(para)
%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2023-10-21
%% This function is to reproduce the SCREENES algorithm.
%%
tic

if ~isfield(para,'rel_tol')
	rel_tol=1e-5;
else
	rel_tol = para.rel_tol;
end 

if ~isfield(para,'lambda')
	lambda = 0.001;
else
	lambda = para.lambda;
end 

FidNus = para.FidNus;
mask = para.mask;

[N1,N2] = size(mask);
FidNus = FidNus(:);

InArg.AM = Fu_downsample(mask,N1,N2);
InArg.AMh = (InArg.AM)';
InArg.y = [real(FidNus);imag(FidNus)];
InArg.ym = length(InArg.y);
InArg.lambda = lambda;
InArg.rel_tol = rel_tol;
InArg.xn = N1*N2*2;

ymax = max(InArg.y);
InArg.y = InArg.y./ymax;

x = tnipm(InArg);
x = x * ymax;
xc = x(1:length(x)/2)+1i*x(length(x)/2+1:end);
final_spectrum=reshape(xc,N1,N2);
Rec_X = final_spectrum;
Rec_X = real(Rec_X);
Rec_X = max(Rec_X,0);
Rec_X = Rec_X./max(max(Rec_X));
Rec_X = fftshift(fftshift(Rec_X,1),2);
rec_time = toc 
fprintf('Completing, time consumes %0.2f sec \n',rec_time);
end
