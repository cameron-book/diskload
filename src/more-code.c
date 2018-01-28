  phase_ind = b>pi/2 | b<0;
  mone_ind= m==1;

// Periodicity for n>1
//
// This is not documented anywhere, but plotting the function in Mathematica
// shows clearly that it's (anti-)symmetric around zero and pi/2.

  phasen_ind=phase_ind & n>1;
  if any(phasen_ind)
  
  mm = m(phasen_ind);
  bb = b(phasen_ind);
  nn = n(phasen_ind);
  
  if any(b>pi/2 & b<pi)
    
    cc=bb-pi/2;
    P(phasen_ind)=conj(-elliptic3x(pi/2-cc,mm,nn));
    
  elseif any(b>pi)
    
    P(phasen_ind)=conj(elliptic3x(bb-pi,mm,nn));
    
  elseif any(b<-pi/2 & b>-pi)
    
    cc=-pi/2-bb;
    P(phasen_ind)=conj(-elliptic3x(-pi/2+cc,mm,nn));
    
  elseif any(b<pi)
    
    P(phasen_ind)=conj(elliptic3x(bb+pi,mm,nn));
    
  end


% Periodicity for n<1:
% http://functions.wolfram.com/EllipticIntegrals/EllipticPi3/04/02/03/0001/
ind = phase_ind & ~mone_ind & n<1;
if any(ind)
  
  mm = m(ind);
  bb = b(ind);
  nn = n(ind);
  
  phi = mod(bb+pi/2,pi)-pi/2;
  a = round(bb./pi);
  P(ind) = 2.*a.*elliptic3c(mm,nn) + sign(phi).*elliptic3x(abs(phi),mm,nn);
  
end


M=(1./sin(b)).^2; %critical value which goes from real inputs to complex inputs
if any(m>M)
  warning('elliptic123:PI_bmn_large_m','Cannot calculate elliptic3(b,m,n) with m greater than the critical value.')
end

% Special case m==n: http://dlmf.nist.gov/19.6#E13
mnequal_ind = m==n & ~phase_ind;
if any(mnequal_ind)
  
  bb=b(mnequal_ind);
  mm=m(mnequal_ind);
  nn=n(mnequal_ind);
  
  [FF,EE]= elliptic12x(bb,mm);
  
  P(mnequal_ind)=(1./(1-mm)).*(EE-((mm./sqrt(1-mm.*(sin(bb)).^2)).*sin(bb).*cos(bb)));
  
end

% Reciprocal-modulus transformation: http://dlmf.nist.gov/19.7#E4
mgnl_ind = m>1 & n<1 & ~phase_ind;
if any(mgnl_ind)
  
  bb=b(mgnl_ind);
  mm=m(mgnl_ind);
  nn=n(mgnl_ind);
  
  P(mgnl_ind)=1./sqrt(mm).*elliptic3ic(asin(sqrt(mm).*sin(bb)),1./mm,nn./mm);
  
end

% Imaginary-modulus transformation: http://dlmf.nist.gov/19.7#E5
mlnl_ind = m~=n & m<0 & n<1 & ~phase_ind;
if any(mlnl_ind)
  
  bb=b(mlnl_ind);
  mm=m(mlnl_ind);
  nn=n(mlnl_ind);
  
  t=asin((sin(bb).*sqrt(1-mm))./sqrt(1-mm.*(sin(bb)).^2));
  
  P(mlnl_ind)=1./((nn-mm).*sqrt(1-mm)).*(-mm.*elliptic12x(t,-mm./(1-mm))+nn.*elliptic3ic(t,-mm./(1-mm),(nn-mm)./(1-mm)));
  
end

% Normal ranges:
mnormnl_ind=n<1 & ~phase_ind & m>=0 & m<=1;
if any(mnormnl_ind)
  
  bb=b(mnormnl_ind);
  mm=m(mnormnl_ind);
  nn=n(mnormnl_ind);
  
  P(mnormnl_ind)=elliptic3ic(bb,mm,nn);
  
end

% Refer to 17.7.8 in Abramowitz:
ng_ind=n>1 & m<n & ~phase_ind; %case where n>1 but m<n
if any(ng_ind)
  
  bb=b(ng_ind);
  mm=m(ng_ind);
  nn=n(ng_ind);
  
  N=mm./nn;
  P1=sqrt(((nn-1).*(1-mm./nn)));
  D=sqrt(1-mm.*(sin(bb)).^2);
  
  P(ng_ind)=-elliptic3x(bb,mm,N)+elliptic12x(bb,mm)+(1./(2.*P1)).*log((D+P1.*tan(bb)).*(D-P1.*tan(bb)).^-1);
  
end

if any(n>1 & ~phase_ind & m>n)
  warning('elliptic123:PI_bmn_large','Cannot calculate elliptic3(b,m,n) for n>1 and m>n.')
end

end
  
}

d[n_,m_] := EllipticPi[n, m] - EllipticPi[n/m,ArcSin[Sqrt[m]],1/m]/Sqrt[m]

d[n_,m_] := EllipticPi[n, m] - EllipticPi[n/m,ArcSin[Sqrt[m]],1/m]/Sqrt[m]

d[n_,m_] := EllipticPi[n, m] - EllipticPi[n/m,ArcSin[Sqrt[m]],1/m]/Sqrt[m]                        

