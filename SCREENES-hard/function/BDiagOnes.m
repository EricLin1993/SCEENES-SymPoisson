function [ B ] = BDiagOnes( N1,len )
%% Author: Enping Lin
%% Contact: enpinglin@qq.com
%% Date: 2023-09-30
%% This script is to find diagonal peaks area.
%%
B=ones(N1);
U=triu(ones(N1),len/2);
L=tril(ones(N1),-len/2);
B=B-U;
B=B-L;
end