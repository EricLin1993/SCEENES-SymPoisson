function res = mtimes(a,b)

if isa(a,'FuSym_downsample_soft') == 0
    error('In  A.*B only A can be Fu_downsample_soft');
end

if a.adjoint
    bA = b(1:length(b)/4) + 1i*b(length(b)/4+1:length(b)/2); 
    bT = b(length(b)/2+1:length(b)*3/4) + 1i*b(length(b)*3/4+1:end); 
%     bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end); 
%     fid2d = zeros(size(a.mask));
%     fid2d(logical(a.mask)) = bc;
    bA = reshape(bA,size(a.mask));
    bT = reshape(bT,size(a.mask));
    fid2d = bA;
    fid2d(a.mask==0) = 0;
    [r,c] = size(fid2d);
    fid2d = fft2(fid2d)/sqrt(r*c);
%     fidvec = Sym2DToVecPer(fid2d);
    fidvec = fid2d(:);
    
    Tbc = Trans_op(bT);
    Tbc = Tbc*sqrt(a.lambda2);
    Tbc = Tbc(:);
    resvec = fidvec + Tbc;
    res = [real(resvec);imag(resvec)];  
else
%     n = size(a.mask,2); 
    bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end);
%     bc = Vec2Sym2DPer(bc,n);
    bc = reshape(bc,size(a.mask));
    [r,c] = size(bc);
    Spec2d = ifft2(bc)*sqrt(r*c); % FX
%     Spec2d_ds = Spec2d(logical(a.mask));% PFX
    Spec2d_ds = Spec2d;
    Spec2d_ds(a.mask==0) = 0;% PFX
    Spec2d_ds = Spec2d_ds(:);
    
    Tbc = Trans_op(bc);
    Tbc = Tbc*sqrt(a.lambda2);
    Tbc = Tbc(:);
    res = [real(Spec2d_ds);imag(Spec2d_ds);real(Tbc);imag(Tbc)];
end