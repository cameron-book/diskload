function LNx = extrapELLNs(LN,Nmax)
%extrapLNs extrapolate the load Love numbers via asymptotic approximation
% This function assumes the the elastic load Love numbers h,k,l complete 
% through degree nmax were computed directly, and that the asymptotic 
% limit has been reached, i.e. for any n > nmax, h(n)=h(nmax), 
% n k(n)= nmax k(nmax), and n l(n)= nmax l(nmax). The function then
% computes h,k and l thru degree Nmax.  Note that this function assumes
% that vectors h,k,l begin with h0,k0,l0 and that length(h)=length(k)=
% length(l)=nmax+1. So, the user must add the zero degree LNs if they
% are not already present before calling this function.
%
% USAGE:   LNx = extrapELLNs(LN,Nmax)
%
% INPUT:
%    LN   the eastic loading Love numbers (LN), a structure with fields
%           LN.h  (nmax+1)-vector containing LN h for degrees 0:nmax
%           LN.k  (nmax+1)-vector containing LN k for degrees 0:nmax
%           LN.k  (nmax+1)-vector containing LN l for degrees 0:nmax
%  Nmax   the desired (larger) value of nmax. Nmax is a scalar.
%
% OUTPUT:
%   LNx   a structure containing the extended LN vectors, with fields
%           LNx.h  (Nmax+1)-vector containing LN h for degrees 0:Nmax
%           LNx.k  (Nmax+1)-vector containing LN k for degrees 0:Nmax
%           LNx.k  (Nmax+1)-vector containing LN l for degrees 0:Nmax


% version 1.0              Michael Bevis            30 March 2017

[n1,n2]=size(LN.h); 
if min([n1 n2])~=1, error('LN.h must be a vector'); end
if n1>n2, icolvec=1; else icolvec=0; end
h=LN.h'; k=LN.k'; l=LN.l';
% check for illegal conditions
m=length(h);
if length(k)~=m || length(l)~=m
    error('Love number vectors do not have the same length')
end

% check for 0-order LNs
if( (l(1)*k(1)) ~= 0 )
    error('n=0 Love numbers l_0 and k_0 must be zero');
end

% check Nmax
if numel(Nmax)~=1 | rem(Nmax,1)~=0
    error('input argument nmax should be scalar (integer)')
end

nn=length(h);
nmax=nn-1;
n=0:nmax;

if Nmax<=nmax
    error('Nmax has to be larger than nmax')
end

h0=h(1);  k0=k(1); l0=l(1);
h(1)=[];  k(1)=[]; l(1)=[];

N=Nmax;
nplus=N-nmax;
hplus=h(end)*ones(1,nplus);
kplus=nmax*k(end)./[nmax+1:N];
lplus=nmax*l(end)./[nmax+1:N];

LNx.h=[h0 h hplus];
LNx.k=[k0 k kplus];
LNx.l=[l0 l lplus];

% if LNs arrived as column (row) vectors, return them the same way
if icolvec==1
    LNx.h=LNx.h(:);
    LNx.k=LNx.k(:);
    LNx.l=LNx.l(:);
end
