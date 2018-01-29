clear 

% Load LNs used in Bevis et al. (2016) GJI
fid = fopen( 'REF_6371_loading_love_numbers_0_40000.txt', 'r' );
data = textscan( fid, '%d %f %f %f', 'headerlines', 1 );
LN.h = data{2};
LN.l = data{3};
LN.k = data{4};

% note: the vectors h,k,l contain the elastic loading Love numbers for
% degrees zero to n       (h0,l0,k0 must be included)

% Set some parameters and constants

% specify angular distances of stations from the disk center
t1 = logspace(-4,log10(175),950);
t2=180-fliplr(logspace(-4,log10(5),51));
theta=union(t1,t2);
% now describe the load   
alpha = 0.1;                           % Disk radius (degrees)
Tw    = 1;                             % Disk height (equivalent water height, m)
icomp = 0;                             % choose imass,0 or 1 (uncompensated/compensated load)
%
if icomp==1
    fprintf('\ninvoking a globally compensated load (icomp=1)\n')
else
    fprintf('\ninvoking an uncompensated load (icomp=0)\n')
end

EM.Re=6371;            % Radius of the Earth (km)
EM.grav=9.8046961;     % Surface gravity (m/s2)

%% Solve the problem using three different values of nmax (4K, 40K, 400K)

fprintf(1,'       NMAX       COMPUTE TIME\n')
for i=1:3
   if i==1
       nmax=4000;
   elseif i==2
       nmax=40000;
   else
       nmax=400000;
       LN = extrapELLNs(LN,nmax);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   tic;
   [U,V,G]= sphericaldiskload(alpha,icomp,theta,Tw,nmax,LN,EM);
   et=toc;
   fprintf(1,'     %6i  ',nmax);
   fprintf(1,'     %8.3f\n',et);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   eval(['save benchmark',int2str(i),'  alpha icomp theta nmax LN EM', ...
        ' U V G ;'])
end

 %%  Save the results from the third run in ascii format too
 fid=fopen('Benchmark.txt','w');
 fprintf(fid, ['   theta (deg) '    '         U (mm)      ',...
               '       V (mm)      ' '      G (mm)      ' '\n'] );
 for j=1:length(theta)
     fprintf(fid,'%15.11f  ',theta(j));
     fprintf(fid,'%17.10e  ',[U(j) V(j) G(j)]);
     fprintf(fid,'\n');
 end
 fclose(fid);
 
