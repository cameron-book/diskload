function Pn = legendre_step(n,x,Pnm1,Pnm2)
%legendre_step  evaluate Pn(x) from Pn-1(x) and Pn-2(x) via a recursion
% relation , where x can be scalar or vector, but x is restricted to 
% the range [-1,1]. 
%
% USAGE:       
%              Pn = legendre_step(n,x,Pnm1,Pnm2)
%
% INPUT:   
%    Pnm1     the Legendre polynomial Pn-1(x) for scalar or vector x
%    Pnm2     the Legendre polynomial Pn-2(x) for scalar or vector x.
%             Must have the same size as Pnm1.
%
% OUTPUT:
%     Pn      the Legendre polynomial Pn(x) 
%
% Normally one will begin the recursion by setting Pnm1 = P1(x) = x and 
% Pnm2 = P0(x) = ones(size(x)), and then using this function to 
% computing Pn = P2(x). Having computed P2 from P1 and P0, it is possible
% to call this function to compute P3 from P2 and P1, then to compute 
% P4 from P3 and P2, and so on.
%
% This routine is useful when working with very large vector arguments x
% and large values of the maximum degree nmax, since unlike legs.m this
% function does not store in memory the results for all previous steps 
% of the recursion relation.
%
% This function was developed for use by function sphericaldiskload.m
%
% Dependencies:  None

% version 1.0               Michael Bevis              21 March 2017 
if size(Pnm1)~=size(Pnm2) 
    error('input arguments Pnm1 and Pnm2 do not have the same size')
end
rn=1./n;
Pn= (2.0-rn).*x.*Pnm1 - (1.0-rn)*Pnm2;

% suppress numerical instabilities that occur when x= +1 or x= -1 exactly
k=find(x==1 | x==-1);
if ~isempty(k)
    Pn(k)=sign(Pn(k));
end