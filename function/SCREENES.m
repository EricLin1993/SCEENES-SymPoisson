function PeakSave = SCREENES(para)
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
InArg.xn = para.xn;

[N1,N2] = size(mask);
FidNus = FidNus(:);

InArg.AM = FuSym_downsample(mask);
InArg.AMh = (InArg.AM)';
InArg.y = [real(FidNus);imag(FidNus)];
InArg.ym = length(InArg.y);
InArg.lambda = lambda;
InArg.rel_tol = rel_tol;

PeakSave = zeros(N1,N2);
itn = 2;

for iter = 1:itn
    InArg.y = [real(FidNus);imag(FidNus)];
    ymax = max(InArg.y);
    InArg.y = InArg.y./ymax;
    % profile on
    x = l1_ls(InArg);
    % profile viewer
    x = x*ymax;
    xc = x(1:length(x)/2)+1i*x(length(x)/2+1:end);
    [ X_Symetric ] = Vec2Sym2DPer( xc,N1 );
    Temp = abs(X_Symetric);
    Temp = Temp./max(Temp(:));
    PSaveF = Temp>0.1; % peak save
   
    if iter == itn
      PeakSave = PeakSave+ X_Symetric;
      break;
    else    
      PeakSave(PSaveF) = PeakSave(PSaveF) + X_Symetric(PSaveF);
    end
    
    PeakSaveFid = ifft2(PeakSave)*sqrt(N1*N2);
    PeakSaveFid(mask==0)= [];
    FidNus = FidNus - PeakSaveFid(:);
end

PeakSave = real(PeakSave);
PeakSave = max(PeakSave,0);
PeakSave = PeakSave./max(max(PeakSave));
rec_time = toc;
fprintf('Completing, time consumes %0.2f sec \n',rec_time)
end
