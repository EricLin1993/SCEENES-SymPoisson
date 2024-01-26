%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2024-01-26
%% This script is to reproduce 15N-15N symmetrical spectra reconstruction using SCREENES-hard and SCREENES-soft.
%% The spectra data are obtained from the 620×620 fully sampled 15N-15N plane projected data matrix of 
%% the 3D (H)N(COCO)NH spectrum of a-synuclein.
%%
clear classes
clear,close all,clc
addpath('data')
addpath('function')
addpath('mask')
addpath('operator')
load('Ny_Nz.mat')

%% fully sampled data
[N1,N2] = size(Spec2D_Ideal);
ppm1 = linspace(107.8,131.6,N1);
ppm2 = linspace(107.8,131.6,N2);
figure,contour(ppm1,ppm2,(Spec2D_Ideal),(linspace(0.005,1,100)))
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('Full Sampling')

%% mask generation
load mask_Ny_Nz.mat

%% NUS data generation
soft = 0;
Fid = ifft2(Spec2D_Ideal);
FidNus = Fid;

if soft
    FidNus(mask==0) = 0;
    FidNus = FidNus(:); % soft
else
    FidNus(mask==0) = [];
end

%% SCREENES
rel_tol = 1e-5;
temp = cumsum(1:N1);
if soft
    xn = 2*length(mask(:)); % soft
else
    xn = 2*temp(end);
end

lambda = 0.001; % 0.001
lambda2 = 0.003; % 0.003
para.FidNus = FidNus;
para.mask = mask;
para.rel_tol = rel_tol;
para.xn = xn;

para.lambda = lambda;
    
if soft
    para.lambda2 = lambda2;
    Rec_X_Symetric = SCREENES_soft(para); % soft
else
    Rec_X_Symetric = SCREENES(para);
end

figure,contour(ppm1,ppm2,Rec_X_Symetric,linspace(0.005,1,100));
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('SCREENES Reconstruction')

%% diagonal peak analysis
RecS = real(Rec_X_Symetric)./max(max(real(Rec_X_Symetric))); % normalization
IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
[pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.07); % find diagonal peaks
bnd1 = floor(locs-w/2);
bnd2 = ceil(locs+w/2); % find diagonal peak range

for it =1:length(locs) % compute peak intensity
    IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
    RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
end

RLNE_Diag = norm(IdeS_diagPInt-RecS_diagPInt)/norm(IdeS_diagPInt);

%% cross peak analysis
IdePeak = [];
RecPeak = [];
[ B ] = BDiagOnes( N1,12 ); % find diagonal peaks area
DiagFlag = logical(  B  );
IdeS_RevDiag = IdeS.*(~DiagFlag);
RecS_RevDiag = RecS.*(~DiagFlag);
bwm = ((IdeS_RevDiag+IdeS_RevDiag.')/2>0.005);
bwm = bwareaopen(bwm,1);
[L,num] = bwlabel(bwm);

for it = 1: num
    Temp = find(L==it);
    IdePeak(it) = sum(IdeS_RevDiag(Temp ),'all');
    RecPeak(it) = sum(RecS_RevDiag(Temp ),'all');
end    

RLNE_Cross = norm(IdePeak-RecPeak)/norm(IdePeak);

%% loc_Ideal
[temp_IIdeal,temp_locIdeal] = findpeaks(IdeS(:),'minpeakheight',0.005);
pIdeal = 0*IdeS;
pIdeal(temp_locIdeal) = 1;
[temp_rIdeal,temp_cIdeal] = find(pIdeal);
rIdeal = ppm2(temp_rIdeal)';
cIdeal = ppm1(temp_cIdeal)';
locIdeal = [rIdeal cIdeal];

%% loc_Rec
[temp_IRec,temp_locRec] = findpeaks(RecS(:),'minpeakheight',0.005);
pRec = 0*RecS;
pRec(temp_locRec) = 1;
[temp_rRec,temp_cRec] = find(pRec);
rRec = ppm2(temp_rRec)';
cRec = ppm1(temp_cRec)';
locSCREENES = [rRec cRec];

Rec_SCREENES = RecS;
