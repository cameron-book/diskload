function [U,V,G]=sphericaldiskload(alpha,icomp,theta,w,nmax,LN,EM)
%DSPHERICALDISKLOAD geoelastic response to a uniform spherical disk load
% This function computes the response to a uniform surface pressure load 
% imposed in a disc of angular radius alpha. The elastic response is found 
% at one or more stations located on the surface of the earth at angular 
% distance(s) theta from the center of the disc load.
% The elastic response is computed using user-supplied elastic loading 
% Love numbers (h,k,l) generated using a specific elastic structure model
% for the earth. If three output arguments are invoked, this function
% also computes the change in the height of the geoid at each station.
% The pressure load imposed within the disk is expressed in terms of the 
% equivalent depth height, thickness) of liquid water (density=1000 kg/m3).
%
% USAGE:
%   [u,v] = sphericaldiskload(alpha,icomp,theta,w,nmax,LN);
%   [u,v] = sphericaldiskload(alpha,icomp,theta,w,nmax,LN,EM);
% [u,v,g] = sphericaldiskload(alpha,icomp,theta,w,nmax,LN);
% [u,v,g] = sphericaldiskload(alpha,icomp,theta,w,nmax,LN,EM);
%
% INPUT VARIABLES:
%   alpha   angular radius of disk stated in degrees (scalar)
%   icomp   switch for a compensated (1) / uncompensated (0) disc
%               load (scalar, integer)
%  theta   angular distances of stations from disc center in degrees
%      w   pressure imposed on the surface within the spherical disk
%             expressed as the height or depth of equivalent water load (m)
%   nmax   maximum harmonic degrees of the expansion to be used. nmax
%          may be a scalar or a vector with multiple truncation points.
%          If nmax=[], then nmax will be set to 5*(360/alpha) if the LNs
%          are complete thru that degree,so as to satisfy Bevis et al. 
%          (2016)s Rule of Thump (ROT) for a good approximation. Otherwise
%          nmax will be set to the highest value of n in LN.
%    LN   the elastic loading Love numbers (LN), a structure with fields
%           LN.h  (n+1)-vector containing Love number h for degrees 0:n
%           LN.k  (n+1)-vector containing Love number k for degrees 0:n
%           LN.l  (n+1)-vector containing Love number l for degrees 0:n
%    EM   (optional) spherical Earth model structure with fields
%             EM.r    radius of the Earth in km. Defaults to 6371.0 km
%            EM.grav  surface gravity in m/s2.   Defaults to 9.8046961 m/s2
%          If EM is input, both fields must be assigned values.
%
% OUTPUT VARIABLES:
%     u     radial or 'vertical' elastic displacement (mm)
%     v     tangential or 'horizontal' elastic displacement (mm)
%     g     geoid change (mm)
%     
%     Size of output vectors is [ length(nmax) length(theta) ]
%     e.g.,   u(i,:) represents U vs theta at the i-th value of nmax
%             u(:,j) represents U vs nmax at the j-th value of theta
%
% NOTES:
% (1) All elements of nmax must be <= n, the maximum order provided for
%     the elastic loading Love numbers
% (2) input w can be positive or negative allowing the user to model
%     the incremental response to incremental pressure changes
% (3) It is easy to switch from the state to the rate problem. If input
%     h is actually the rate of change of the load (in m/yr w.e.), then 
%     outputs u,v and g will become velocities (in mm/yr).
%  
% REFERENCE:
% This function is a modified verison of function diskload.m associated 
% with the publication 
%    Bevis, M., Melini, D., Spada, G., 2016. On computing the geoelastic 
%    response to a disk load, Geophys. J. Int. 205 (3), 1,804-1,812,
%    doi:10.1093/gji/ggw115
%
% This function is an 'out of core' version of diskload.m, meaning that it
% makes far less use of memory in the event of a large value of nmax,
% because it does not precompute ans store all the Legendre polynomials/
% The argument nmin used by function diskload.m does not appear in this
% function which assumes that nmin = 0.
%
% DEPENDENCIES:
% This function calls functions legs.m and legendre_step.m
%
 
%  version 1.0             Michael Bevis            28 March 2017

% Define or compute some constants

ggg=6.67384e-11;        % Newton's constant (SI units) 
if nargin<7
    Re=6371;            % Radius of the Earth (km)
    grav=9.8046961;     % Surface gravity (m/s2)
else
    Re=EM.Re;
    grav=EM.grav;
end
radiusm=Re*1e3;     % Radius of the Earth (m)
rhow = 1000;            % density of pure water(kg/m^3) 
rhoear=3.0*grav/4.0/ggg/pi/radiusm;       % Average Earth density (kg/m^3)
from_m_to_mm = 1000;

h=LN.h;
k=LN.k;
l=LN.l;

% check for illegal conditions
m=length(h);
if length(k)~=m || length(l)~=m
    error('Love number vectors do not have same length')
end
if ~isempty(nmax) & nmax>(m-1)
    error('nmax exceeds the lengths of the Love Number vectors')
end

% check for 0-order LNs
if( (l(1)*k(1)) ~= 0 )
    error('n=0 Love numbers l_0 and k_0 must be zero');
end

% self-select nmax is this is desired
if isempty(nmax)
    nmx=ceil(5*(360/alpha));
    if nmx > (m-1)
        nmax=(m-1);
    else
        nmax=nmx;
    end
end

% Convert theta to a row vector
dim = size(theta);
if ( dim(2)==1 ), theta = theta'; end

% Computing the harmonic coefficients of the load 
% Vectors are "offset-indexed", i.e.
% P_n = leg(n+1), sigma_n = sigma(n+1)

n_save=nmax;         % the original nmax vector
nmax=nmax(end);
nn=length(n_save);   % the number of different nmax values invoked

calpha=cosd(alpha);
leg   = legs( nmax+1, calpha );   % compute all Pn(calpha) for n=0:nmax
sigma = nan([nmax+1 1]);          % preallocate space

switch icomp
    case 0        % Uncompensated disc load, eq. (7)
        for n=0:nmax
            if( n==0 ), sigma(n+1) = 0.5*(1-calpha); end
            if( n>0  ), sigma(n+1) = 0.5*(-leg(n+2)+leg(n)); end
        end
    case 1        % Compensated disc load, eq. (8)
        for n=0:nmax
            if( n==0 ), sigma(n+1) = 0; end
            if( n>0  ), sigma(n+1) = -(leg(n+2)-leg(n)) ... 
                                                / (1+calpha); end
        end
end

nt=length(theta);
u = zeros(1,nt);
v =  zeros(1,nt);
U=zeros(nn,nt);
V=zeros(nn,nt);
if nargout>2
    g =  zeros(1,nt);
    G =  zeros(nn,nt);
end


x   = cosd(theta);
s   = sind(theta);
idx = (abs(x) == 1);

% Extract the n=0 terms from h,k,l
h0=h(1); k0=k(1); l0=l(1);
h(1)=[]; k(1)=[]; l(1)=[];

% Extract the n=0 terms from sigma
sigma0=sigma(1);
sigma(1)=[];

% Compute the n=0 terms
P0=ones(size(theta));
u = u + h0 * sigma0;
if nargout>2
     g = g + sigma0;
end

% Compute the n=1 terms
P1=x;
u = u + h(1) * sigma(1) * x /3;
dP1=(-1./s).*(P0 - x.*P1);
dP1(idx) = 0.;
v = v + l(1) * sigma(1) * dP1/3;
if nargout>2
     g = g + (1+ k(1)) * sigma(1) * x /3;
end

Pnm2=P0;
Pnm1=P1;
for n=2:nmax
    Pn=legendre_step(n,x,Pnm1,Pnm2);
    dPn=(-n./s).*(Pnm1 - x.*Pn);
    dPn(idx) = 0.;
        
    ampl = sigma(n)/(2*n+1);
    u  =   u +     h(n)  * ampl * Pn;
    v  =   v +     l(n)  * ampl * dPn;
    if nargout>2
          g  =   g +  (1+k(n)) * ampl * Pn;
    end
    
    % save intermediate results if there are multiple values for nmax
   [lia,ind]=ismember(n,n_save);
   if lia
        U(ind,:)=u;
        V(ind,:)=v;
        if nargout>2
            G(ind,:)=g;
        end
    end
    
    Pnm2=Pnm1;
    Pnm1=Pn;
end

U = U * (3*rhow/rhoear) * w * from_m_to_mm;
V = V * (3*rhow/rhoear) * w * from_m_to_mm;
if nargout>2
    G = G * (3*rhow/rhoear) * w * from_m_to_mm;
end


end

    

        
        

