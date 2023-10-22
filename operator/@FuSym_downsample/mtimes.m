function res = mtimes(a,b)

if isa(a,'FuSym_downsample') == 0
    error('In  A.*B only A can be Fu_downsample');
end

if a.adjoint
    bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end); 
    fid2d=zeros(size(a.mask));
    fid2d(logical(a.mask))=bc;
    [r,c]=size(fid2d);
    fid2d=fft2(fid2d)/sqrt(r*c);
    fidvec = Sym2DToVecPer(fid2d);
    fidvec = fidvec(:);
    res = [real(fidvec);imag(fidvec)];  
else
    n = size(a.mask,2); 
    bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end);
    bc = Vec2Sym2DPer(bc,n);
    [r,c]=size(bc);
    Spec2d=ifft2(bc)*sqrt(r*c);
    Spec2d_ds=Spec2d(logical(a.mask));
    Spec2d_ds = Spec2d_ds(:);
    res = [real(Spec2d_ds);imag(Spec2d_ds)];
end