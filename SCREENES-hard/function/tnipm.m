function [x] = tnipm(InArg)
%% Author: Enping Lin, Ze Fang, Bo Chen
%% Contact: enpinglin@qq.com
%% Date: 2023-10-21
%% This function is to solve the l1-regularized least squares problems of the following form:
%% minimize ||A*x-y||^2 + lambda*sum|x_i|.
%% The code cites the article "An interior-point method for large-scale
%% l1-regularized least squares" (S. Kim, K. Koh, M. Lustig, S. Boyd, D. Gorinevsky)

%% parameter setting
mu = 2; % updating parameter of t
max_nt_iter = 400; % maximum IPM (Newton) iteration

alpha = 0.01; % minimum fraction of decrease in the objective
beta = 0.5; % stepsize decrease factor
max_ls_iter = 100; % maximum backtracking line search iteration

if ~isfield(InArg,'rel_tol')
	rel_tol=1e-5;
else
	rel_tol = InArg.rel_tol;
end 

if ~isfield(InArg,'lambda')
	lambda = 0.001;
else
	lambda = InArg.lambda;
end 

A = InArg.AM;
At = InArg.AMh;
m = InArg.ym;
n = InArg.xn;
y = InArg.y;

eta = 1e-3;
pcgmaxi = 5000;
x = zeros(n,1);
u = ones(n,1);
t = min(max(1,1/lambda),2*n/1e-3);

f = [x-u;-x-u];

%% result and history variables
pobjs = []; 
dobjs = []; 
sts = [];  
pitrs = []; 
pflgs = [];
dobj = -Inf; 
s = Inf;
pitr = 0 ;
pflg = 0 ;

dxu =  zeros(2*n,1);
diagxtx = 2*ones(n,1);

disp(sprintf('\nSolving a problem of size (m=%d, n=%d), with lambda=%.5e',m,n,lambda)); 
disp('-----------------------------------------------------------------------------');
disp(sprintf('%5s %9s %15s %15s %13s %11s','iter','gap','primobj','dualobj','step len','pcg iters')); 

%% main loop
for ntiter = 0:max_nt_iter
    z = A*x-y;
    
    %% duality gap calculating
    nu = 2*z;
    maxAnu = norm(At*nu,inf);
    
    if (maxAnu > lambda)
        nu = nu*lambda/maxAnu;
    end
    
    pobj  =  z'*z+lambda*norm(x,1);
    dobj  =  max(-0.25*nu'*nu-nu'*y,dobj);
    gap   =  pobj - dobj;

    pobjs = [pobjs pobj]; 
    dobjs = [dobjs dobj]; 
    sts = [sts s];
    pflgs = [pflgs pflg]; 
    pitrs = [pitrs pitr];

    disp(sprintf('%4d %12.2e %15.5e %15.5e %11.1e %8d',ntiter, gap, pobj, dobj, s, pitr)); 

    if (gap/dobj < rel_tol) 
        disp('Absolute tolerance reached.'); 
        return;
    end   
    
    if (s >= 0.5)
        t = max(min(2*n*mu/gap, mu*t), t);
    end

    %% newton step calculating
    q1 = 1./(u+x);
    q2 = 1./(u-x);
    d1 = (q1.^2+q2.^2)/t;
    d2 = (q1.^2-q2.^2)/t;

    %% gradient calculating
    gradphi = [At*(z*2)-(q1-q2)/t; lambda*ones(n,1)-(q1+q2)/t];
    
    %% vectors to be used in the preconditioner calculating
    prb = diagxtx+d1;
    prs = prb.*d1-(d2.^2);

    %% pcg tolerance setting
    normg = norm(gradphi);
    pcgtol = min(1e-1,eta*gap/min(1,normg));
    
    if (ntiter ~= 0 && pitr == 0) 
        pcgtol = pcgtol*0.1; 
    end

    [dxu,pflg,~,pitr,~] = ...
        pcg(@AXfunc_l1_ls,-gradphi,pcgtol,pcgmaxi,@Mfunc_l1_ls,...
            [],dxu,A,At,d1,d2,d1./prs,d2./prs,prb./prs);

    if (pflg == 1) 
        pitr = pcgmaxi; 
    end
    
    dx = dxu(1:n);
    du = dxu(n+1:end);
    
    %% backtracking line search algorithm
    phi = z'*z+lambda*sum(u)-sum(log(-f))/t;
    s = 1.0;
    gdx = gradphi'*dxu;
    for lsiter = 1:max_ls_iter
        newx = x+s*dx;
        newu = u+s*du;
        newf = [newx-newu;-newx-newu];
        if (max(newf) < 0)
            newz = A*newx-y;
            newphi = newz'*newz+lambda*sum(newu)-sum(log(-newf))/t;
            if (newphi-phi <= alpha*s*gdx)
                break;
            end
        end
        s = beta*s;
    end
    if (lsiter == max_ls_iter) 
        break; 
    end
        
    x = newx;
    u = newu;
    f = newf;
end
return;

function [y] = AXfunc_l1_ls(x,A,At,d1,d2,p1,p2,p3)
%% compute AX (PCG)
    n  = length(x)/2;
    x1 = x(1:n);
    x2 = x(n+1:end);

y = [(At*((A*x1)*2))+d1.*x1+d2.*x2; d2.*x1+d1.*x2];% y = hessphi*[x1;x2],
                                                   % where hessphi = [A'*A*2+D1 , D2;
                                                   %                  D2        , D1];

function [y] = Mfunc_l1_ls(x,A,At,d1,d2,p1,p2,p3)
%% compute P^{-1}X (PCG)
    n  = length(x)/2;
    x1 = x(1:n);
    x2 = x(n+1:end);

y = [ p1.*x1-p2.*x2;...
         -p2.*x1+p3.*x2];% y = P^{-1}*x
