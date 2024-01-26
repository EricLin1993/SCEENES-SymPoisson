%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2024-01-26
%% This script is to reproduce statistical analysis result of four sampling schedule
%% (Random, Sym Random, 2D Poisson and Sym 2D Poisson sampling schedules) by SCREEN, using
%% a 15N-15N symmetrical spectra NUS data. 
%% The spectra data are obtained from the 620×620 fully sampled 15N-15N plane projected data matrix of 
%% the 3D (H)N(COCO)NH spectrum of a-synuclein.
%%
clear,close all,clc
addpath('data')
addpath('operator')
addpath('function')
load('Ny_Nz.mat')

%% read data
smodel = {'ran','sym','2Dpoisson','sympoisson'};
lambda = [0.01, 0.0075, 0.005, 0.0025, 0.001];
nusr = [0.25, 0.2, 0.15, 0.1, 0.05];
testnum = 50;
rel_tol = 1e-5;
Rlt_CroRlne = cell(testnum,1);

for it = 1:testnum
	Rlt_CroRlne{it} = zeros(length(lambda),length(nusr),length(smodel));
end
Rlt_DiaRlne = Rlt_CroRlne;

[N1,N2]=size(Spec2D_Ideal);

%% statistical data generation
for itr_test = 1:testnum
	for itr_smodel = 1:length(smodel)
        for itr_lambda = 1:length(lambda)
            for itr_nusr = 1:length(nusr)         
               %% NUS sampling
                [ mask ] = SNUSMask( N1,nusr(itr_nusr),smodel{itr_smodel} ); 
                Fid = ifft2(Spec2D_Ideal);
                FidNus = Fid;
                FidNus(mask==0) = [];
                FidNus = FidNus(:);

               %% SCREEN
                para.FidNus = FidNus;
                para.mask = mask;
                para.rel_tol = rel_tol;
                para.lambda = lambda(itr_lambda);

                Rec_X_Symetric = SCREEN(para);
         
                RecS = real(Rec_X_Symetric)./max(max(real(Rec_X_Symetric)));% normalization
                IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
                
               %% diagonal peak analysis
                [pks,locs,w,p] = findpeaks(diag(RecS),'minpeakheight',0.01);  % find diagonal peaks
                bnd1 = floor(locs-w/2);bnd2 = ceil(locs+w/2); % find diagonal peak range
               
                for it = 1:length(locs) % compute peak intensity
                   IdeS_diagPInt(it) = sum(sum(IdeS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
                   RecS_diagPInt(it) = sum(sum(RecS(bnd1(it):bnd2(it),bnd1(it):bnd2(it))));
                end
                
                RLNE_Diag = norm(IdeS_diagPInt-RecS_diagPInt)/norm(IdeS_diagPInt);
                Rlt_DiaRlne{itr_test}(itr_lambda,itr_nusr,itr_smodel) = RLNE_Diag;
                
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
                
                RLNE_Cross = norm(IdePeak-RecPeak)/norm(IdePeak);
                Rlt_CroRlne{itr_test}(itr_lambda,itr_nusr,itr_smodel) = RLNE_Cross;
            end    
        end    
	end    
end

%% statistic analysis
DiaRlneArray = zeros(length(lambda),length(nusr),length(smodel),testnum);
CroRlneArray = zeros(length(lambda),length(nusr),length(smodel),testnum);

for it = 1:testnum
    DiaRlneArray(:,:,:,it) = Rlt_DiaRlne{it};
	CroRlneArray(:,:,:,it) = Rlt_CroRlne{it};
end  

Mean_DiaRlne = mean(DiaRlneArray,4);
Mean_CroRlne = mean(CroRlneArray,4);
Std_DiaRlne = std(DiaRlneArray,1,4);
Std_CroRlne = std(CroRlneArray,1,4);

save('SCREEN_data.mat','Mean_DiaRlne','Mean_CroRlne','Std_CroRlne','Std_DiaRlne')