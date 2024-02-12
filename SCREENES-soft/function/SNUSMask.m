function [ mask ] = SNUSMask( N,nusr,smodel )
%% Author: Enping Lin
%% Contact: enpinglin@qq.com
%% Date: 2023-09-30
%% This function is to generate the sampling mask.
%%
switch smodel
	case '2Dpoisson'    
        inputArg.dim=2;
        inputArg.n1=N;
        inputArg.n2=N;
        inputArg.n3=1;
        inputArg.sr=nusr;
        inputArg.sineportion = 1;
     
        [mask] = PG_HS_LEP(inputArg);

	case 'sym'   
        c = N;
        n = c*(c+1)/2;
        m = zeros(n,1); 
        l = round(nusr*N*N);
        ind = randperm(n,l);
        m(ind)=1;
        m = logical(m);
        upf = zeros(n,1);
        upf(rand(size(upf))<0.5)=1;
        [ mask ] = SymMask( m,upf,c);
	case 'ran'
        l = round(nusr*N*N);
        mask = zeros(N,N);
        ind = randperm(N*N,l);
        mask(ind)=1;
	case 'sympoisson'    
        c = N;
        n = c*(c+1)/2;
        
        m = zeros(n,1);
        [v,k] = Poisson1D(1,round(N*N*nusr),n);
        m(v)=1;
        m = logical(m);
    
        upf = zeros(n,1);
        upf(rand(size(upf))<0.5)=1;
        [ mask ] = SymMask( m,upf,c);

	otherwise  
        disp('smodel is set ran in default...')
        n = N*N;
        l = round(nusr*n);
        mask = zeros(N,N);
        ind = randperm(n,l);
        mask(ind)=1;
	end
end

