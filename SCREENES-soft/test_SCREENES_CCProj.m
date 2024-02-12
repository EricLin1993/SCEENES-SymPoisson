%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2023-10-21
%% This script is to reproduce 13C-13C symmetrical spectra reconstruction using SCREENES.
%% The spectra data are obtained from 512*512 13C-13C projection of 4D methyl HMQC-NOESY-HMQC spectrum.

%%
clear,close all,clc
addpath('data')
addpath('function')
addpath('mask')
addpath('operator')
load('CCProj_spec.mat')

%% fully sampled data
[N1,N2] = size(Spec2D_Ideal);
ppm1 = linspace(14,27.2,N1);
ppm2 = linspace(14,27.2,N2);
figure,contour(ppm1,ppm2,(Spec2D_Ideal),(linspace(0.005,1,50)))
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('Full Sampling')

%%
smodel = {'ran','sym','2Dpoisson','sympoisson'};
lambda = [0.01, 0.01, 0.01, 0.01];

for iter = 1:length(smodel)
    %% mask generation
    if iter == 1
        load mask_ran.mat
    elseif iter == 2
        load mask_sym.mat
    elseif iter == 3
        load mask_2Dpoisson.mat
    else
        load mask_sympoisson.mat
    end

    %% NUS data generation
    Fid = ifft2(Spec2D_Ideal);
    FidNus = Fid;
    FidNus(mask==0) = [];

    %% SCREENES
    rel_tol = 1e-5;
    temp = cumsum(1:N1);
    xn = 2*temp(end);

    para.FidNus = FidNus;
    para.mask = mask;
    para.rel_tol = rel_tol;
    para.xn = xn;
    para.lambda = lambda(iter);

    Rec_X_Symetric = SCREENES(para);
    
    figure,contour(ppm1,ppm2,Rec_X_Symetric,linspace(0.005,1,50));
    set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
    title(['SCREENES Reconstruction(', smodel{iter},')'])
    
    %% diagonal peak analysis
    RecS = real(Rec_X_Symetric)./max(max(real(Rec_X_Symetric))); % normalization
    IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
    [pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.03); % find diagonal peaks
    bnd1 = floor(locs-w/2);
    bnd2 = ceil(locs+w/2); % find diagonal peak range

    for it = 1:length(locs) % compute peak intensity
       IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
       RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
    end

    sl = max(IdeS_diagPInt)*1.1;
    Dia_lx{iter} = polyfit(IdeS_diagPInt,RecS_diagPInt,1); % curve fitting

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

    Cro_lx{iter} = polyfit(IdePeak,RecPeak,1);
end
