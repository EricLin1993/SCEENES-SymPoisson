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
soft = 1;
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
bwm = bwareaopen(bwm,2);
[L,num] = bwlabel(bwm);

for it = 1: num
    Temp = find(L==it);
    IdePeak(it) = sum(IdeS_RevDiag(Temp ),'all');
    RecPeak(it) = sum(RecS_RevDiag(Temp ),'all');
end    

RLNE_Cross = norm(IdePeak-RecPeak)/norm(IdePeak);

%% =========================== NUS-con ========================
%% loc_Ideal
bwm_Ideal = (IdeS>0.005);
bwm_Ideal = bwareaopen(bwm_Ideal,2) ;
[L_Ideal,num_Ideal] = bwlabel(bwm_Ideal);

locIdeal = zeros(num_Ideal,2);
for it = 1: num_Ideal
    Temp_Ideal = find(L_Ideal==it);
    amp_Ideal = sum(IdeS(Temp_Ideal ),'all');
    [IdeS_X, IdeS_Y] = find(L_Ideal==it);
    area_peak_Ideal = length(IdeS_X);
    loc_temp_Ideal = zeros(area_peak_Ideal,2);
    for aa = 1:area_peak_Ideal
        loc_temp_Ideal(aa,:) = IdeS(IdeS_X(aa), IdeS_Y(aa)).*([ppm2(IdeS_X(aa))',ppm1(IdeS_Y(aa))]);
    end
    locIdeal(it,:) = ([sum(loc_temp_Ideal(:,1)),sum(loc_temp_Ideal(:,2))])/amp_Ideal;
end 

%% loc_Rec
bwm_Rec = (RecS>0.005);
bwm_Rec = bwareaopen(bwm_Rec,2) ;
[L_Rec,num_Rec] = bwlabel(bwm_Rec);

locRec = zeros(num_Rec,2);
for it = 1: num_Rec
    Temp_Rec = find(L_Rec==it);
    amp_Rec = sum(RecS(Temp_Rec ),'all');
    [RecS_X, RecS_Y] = find(L_Rec==it);
    area_peak_Rec = length(RecS_X);
    loc_temp_Rec = zeros(area_peak_Rec,2);
    for bb = 1:area_peak_Rec
        loc_temp_Rec(bb,:) = RecS(RecS_X(bb), RecS_Y(bb)).*([ppm2(RecS_X(bb))',ppm1(RecS_Y(bb))]);
    end
    locRec(it,:) = ([sum(loc_temp_Rec(:,1)),sum(loc_temp_Rec(:,2))])/amp_Rec;
end 

results = double(py.utils.nuscon_metrics(py.numpy.array(IdeS), py.numpy.array(RecS), ...
    py.numpy.array(locIdeal), py.numpy.array(locRec), 1));
results(4) = 1 - results(4);