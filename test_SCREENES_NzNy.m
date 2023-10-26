%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2023-10-21
%% This script is to reproduce 15N-15N symmetrical spectra reconstruction using SCREENES.
%% The spectra data are obtained from 512*512 15N-15N projection of the TROSY-NOESY-TROSY spectrum
%% of the main protease of SARS-CoV-2.
%%
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
figure,contour(ppm1,ppm2,(Spec2D_Ideal),(linspace(0.005,1,50)))
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('Full Sampling')

%% mask generation
load mask_NzNy.mat

%% NUS data generation
Fid = ifft2(Spec2D_Ideal);
FidNus = Fid;
FidNus(mask==0) = [];

%% SCREENES
rel_tol = 1e-5;
temp = cumsum(1:N1);
xn = 2*temp(end);
lambda = [0.001, 0.0025, 0.005, 0.0075, 0.01];

para.FidNus = FidNus;
para.mask = mask;
para.rel_tol = rel_tol;
para.xn = xn;


for iter = 1:length(lambda)
    para.lambda = lambda(iter);
    
    Rec_X_Symetric = SCREENES(para);
    
    lam = num2str(para.lambda);
    figure,contour(ppm1,ppm2,Rec_X_Symetric,linspace(0.005,1,50));
    set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
    title(['SCREENES Reconstruction(lambda=', lam,')'])

    %% diagonal peak analysis
    RecS = real(Rec_X_Symetric)./max(max(real(Rec_X_Symetric))); % normalization
    IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
    [pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.07); % find diagonal peaks
    bnd1 = floor(locs-w/2);bnd2 = ceil(locs+w/2); % find diagonal peak range

    for it =1:length(locs) % compute peak intensity
       IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
       RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
    end

    RLNE_Diag(iter) = norm(IdeS_diagPInt-RecS_diagPInt)/norm(IdeS_diagPInt);

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
       IdePeak(it) = sum(IdeS(Temp ),'all');
       RecPeak(it) = sum(RecS(Temp ),'all');
    end    

    RLNE_Cross(iter) = norm(IdePeak-RecPeak)/norm(IdePeak);
end