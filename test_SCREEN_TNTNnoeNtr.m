%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2023-10-21
%% This script is to reproduce 15N-15N symmetrical spectra reconstruction using SCREEN.
%% The spectra data are obtained from 512*512 15N-15N projection of the (H)N(COCO)NH spectrum of a-synuclein.
%%
clear,close all,clc
addpath('data')
addpath('function')
addpath('mask')
addpath('operator')
load('TNTNnoeNtr.mat')

%% fully sampled data
[N1,N2] = size(Spec2D_Ideal);
ppm1 = linspace(100.9,132.5,N1);
ppm2 = linspace(100.9,132.5,N2);
figure,contour(ppm1,ppm2,(Spec2D_Ideal),(linspace(0.003,1,50)))
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('Full Sampling')

%% mask generation
load mask_TNTNnoeNtr.mat

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
AM = Fu_downsample(mask,N1,N2);
AMh = AM';
Rec_X = SCREEN(para);

figure,contour(ppm1,ppm2,Rec_X,linspace(0.003,1,50));
set (gca,'xdir','reverse','ydir','reverse','box','on','linewidth',1.5,'fontsize',14,'fontname','Calibri','fontweight','bold');
title('SCREEN Reconstruction')

%% diagonal peak analysis
RecS = real(Rec_X)./max(max(real(Rec_X))); % normalization
IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
[pks,locs,w,p]=findpeaks(diag(RecS),'minpeakheight',0.05); % find diagonal peaks
bnd1 = floor(locs-w/2);bnd2 = ceil(locs+w/2); % find diagonal peak range

for it =1:length(locs) % compute peak intensity
   IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
   RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
end

sl = max(IdeS_diagPInt)*1.1;
Dia_lx =polyfit(IdeS_diagPInt,RecS_diagPInt,1); % curve fitting
figure,plot(IdeS_diagPInt,RecS_diagPInt,'o','MarkerFaceColor','b','MarkerEdgeColor','b'), xlim([0,sl]), ylim([0,sl])
hold on, line([0,sl],[0*Dia_lx(1)+Dia_lx(2),sl*Dia_lx(1)+Dia_lx(2)],'LineWidth',2,'color','b')
set(gca,'fontsize',12,'linewidth',1.5,'fontweight','bold','fontname','Calibri',...
    'XColor',[0 0 1],'YColor',[0 0 1],'ZColor',[0 0 1])
title('DiagnolPeaks Intensity')

%% cross peak analysis
IdePeak = [];
RecPeak = [];
[ B ] = BDiagOnes( N1,16 ); % find diagonal peaks area
DiagFlag = logical(  B  );
IdeS_RevDiag = IdeS.*(~DiagFlag);
RecS_RevDiag = RecS.*(~DiagFlag);
bwm = ((IdeS_RevDiag+IdeS_RevDiag.')/2>0.003);
bwm = bwareaopen(bwm,1);
[L,num] = bwlabel(bwm);

for it = 1: num
   Temp = find(L==it);
   IdePeak(it) = sum(IdeS(Temp ),'all');
   RecPeak(it) = sum(RecS(Temp ),'all');
end

Cro_lx = polyfit(IdePeak,RecPeak,1);
lmax = max(IdePeak);
fl = lmax*1.1;
figure,plot(IdePeak,RecPeak,'o','MarkerFaceColor','r','MarkerEdgeColor','r'), xlim([0,fl]), ylim([0,fl])
hold on, line([0,fl],[0*Cro_lx(1)+Cro_lx(2),fl*Cro_lx(1)+Cro_lx(2)],'LineWidth',2,'color','r')
set(gca,'fontsize',12,'linewidth',1.5,'fontweight','bold','fontname','Calibri',...
    'XColor',[1 0 0],'YColor',[1 0 0],'ZColor',[1 0 0])
title('CrossPeaks Intensity')

%% weak peaks analysis
weak_IdePeak_Index = find(IdePeak<0.02);

for it = 1: length(weak_IdePeak_Index)
   weak_IdePeak(it) = IdePeak(weak_IdePeak_Index(it));
   weak_RecPeak(it) = RecPeak(weak_IdePeak_Index(it));
end
Weak_lx = polyfit(weak_IdePeak, weak_RecPeak,1);
wlmax = max(weak_IdePeak);
wl = wlmax*1.1;
figure,plot(weak_IdePeak,weak_RecPeak,'o','MarkerFaceColor','r','MarkerEdgeColor','r'), xlim([0,wl]), ylim([0,wl])
hold on, line([0,wl],[0*Weak_lx(1)+Weak_lx(2),wl*Weak_lx(1)+Weak_lx(2)],'LineWidth',2,'color','r')
set(gca,'fontsize',12,'linewidth',1.5,'fontweight','bold','fontname','Calibri',...
    'XColor',[1 0 0],'YColor',[1 0 0],'ZColor',[1 0 0])
title('WeakPeaks Intensity')