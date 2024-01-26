%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2024-01-26
%% This script is to reproduce 1H-1H symmetrical spectra reconstruction using SCREEN.
%% The spectra data are obtained from 128*128 1H-1H projection of the COSY-TOCSY spectrum of azithromycin.
%%
clear classes
clear,close all,clc
addpath('data')
addpath('function')
addpath('mask')
addpath('operator')
load('azithromycin.mat')

%% fully sampled data
[N1,N2] = size(Spec2D_Ideal);
ppm1 = linspace(0.41,5.41,N1);
ppm2 = linspace(0.41,5.41,N2);
figure,contour(ppm1,ppm2,(Spec2D_Ideal),(linspace(0,1,150)))
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('Full Sampling')

%% mask generation 
load mask_azithromycin.mat

%% NUS data generation
Fid = ifft2(Spec2D_Ideal);
FidNus = Fid;
FidNus(mask==0) = [];
FidNus = FidNus(:);

%% SCREEN
rel_tol = 1e-5;
lambda = 0.001;

para.FidNus = FidNus;
para.mask = mask;
para.rel_tol = rel_tol;
para.lambda = lambda;
    
Rec_X = SCREEN(para);
    
lam = num2str(para.lambda);
figure,contour(ppm1,ppm2,Rec_X,linspace(0,1,150));
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('SCREEN Reconstruction')

%% diagonal peak analysis
RecS = real(Rec_X)./max(max(real(Rec_X))); % normalization
IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
[pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.02); % find diagonal peaks
bnd1 = floor(locs-w/2);
bnd2 = ceil(locs+w/2);  % find diagonal peak range

for it =1:length(locs) % compute peak intensity
    IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
    RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
end

RLNE_Diag = norm(IdeS_diagPInt-RecS_diagPInt)/norm(IdeS_diagPInt);

%% cross peak analysis
IdePeak = [];
RecPeak = [];
[ B ] = BDiagOnes( N1,10 ); %find diagonal peaks area
DiagFlag = logical(  B  );
IdeS_RevDiag = IdeS.*(~DiagFlag);
RecS_RevDiag = RecS.*(~DiagFlag);
bwm = ((IdeS_RevDiag+IdeS_RevDiag.')/2>0.02);
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
locSCREEN = [rRec cRec];

Rec_SCREEN = RecS;
