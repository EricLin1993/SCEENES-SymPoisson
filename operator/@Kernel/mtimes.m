function res = mtimes(a,b)

if isa(a,'Kernel') == 0
    error('In  A.*B only A can be Kernel');
end

if a.adjoint
    bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end); 
   
    fid2d=zeros(size(a.mask));
    fid2d(logical(a.mask))=bc;
    fid2d=fft2_norm(fid2d);
    fidvec = Sym2DPer2Vec(fid2d);
    res = [real(fidvec);imag(fidvec)];
    
else
    bc = b(1:length(b)/2) + 1i*b(length(b)/2+1:end); %
    bc = Vec2Sym2DPer(bc, size(a.mask,1));
    Spec2d=ifft2_norm(bc);
    Spec2d_ds=Spec2d(logical(a.mask));
    Spec2d_ds = Spec2d_ds(:);
    res = [real(Spec2d_ds);imag(Spec2d_ds)];
end