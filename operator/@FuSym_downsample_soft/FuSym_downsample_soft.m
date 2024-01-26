function res =FuSym_downsample_soft(mask,lambda2)

res.adjoint = 0;
res.mask = mask;
res.lambda2 = lambda2;
res = class(res,'FuSym_downsample_soft');
