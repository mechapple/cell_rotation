% Determine Anistropic and isotropic dislocation line energy and
% dislocation nucleation form a crack tip.
%
% This script provides the data for the Anisotropic dislocation line energy
% and crack tip dislocation nucleation ARL technical report. Equation
% numbers given in this script correspond to equation numbers given in the
% report.
%
% Anistropic dislocation orientation prefactor tensor determined from
% equations given by Barnett & Swanger (PhysStat.Sol 1971) and Barnett &
% Asaro (JMPS 1972).
%
% KII isotropic and anisotropic crack tip loading factors for dislocation
% nucleation using in Rice's (JMPS 1994) isotropic dislocation nucleation
% model or Sun & Beltz's (JMPS 1994) anisotropic dislocation nucleation
% model.
%
% Lynn Munday, August 2012
%
%close all hidden,clear all, clc
format short
%------MANUAL INPUT--------------------------------------------------------
%
% RDX stable and unstable stacking energy J/m^2=Pa*m
% From Munday thesis (UMd 2011) and Munday et al. (Phil Mag 2012)
y_us=164e-3; % for (010)[100]
y_sf=101e-3;
% y_us=260e-3; % for (001)[100]
% y_sf=206e-3;
%
% y_us=255e-3; % for (011)[100]
% y_sf=140e-3;
%
% y_us=250e-3; % for (021)[100]
% y_sf=187e-3;
%lattice vectors a,b,c or [100][010][001]
% latVec=[6.4012 0 0;...
% 0 10.2542 0;...
% -1.1243 0 7.5143]*1e-10;

latVec=1e-10*dlmread('cell.mat');
chosen=dlmread('chosen.mat');
chosen;
bvec_n=chosen(1:3); %unit burgers vector / slip direction
lvec=chosen(4:6);

[bvec_n' lvec'];

line_dir=latVec*lvec'; %disocation line direction
bvec=latVec*bvec_n'; %burgers vector

[bvec line_dir];
% RDX properties (MPa) Table 3.6, p79 of Munday thesis (UMd 2011) and
% Munday et al. (JPCB 2010)
% Voigt Stiffness Coefficients from Munday Dissertation
% C11=22.2e3;
% C22=23.9e3;
% C33=23.4e3;
% C44=9.2e3;
% C55=11.1e3;
% C66=10.1e3;
% C23=13.0e3;
% C13=13.2e3;
% C12=9.6e3;
% C15=-0.1e3;
% C25=4.7e3;
% C35=1.6e3;
% C46=2.5e3;
%--------END MANUAL INPUT--------------------------------------------------
% % Voight notation
% Cv=[C11 C12 C13 0 C15 0;...
% C12 C22 C23 0 C25 0;...
% C13 C23 C33 0 C35 0;...
% 0 0 0 C44 0 C46;...
% C15 C25 C35 0 C55 0;...
% 0 0 0 C46 0 C66];

C = 1000.0*dlmread('sewell.mat'); 
Q=dlmread('rotation.mat');
Om = Q

K1 = [Om(1,1)^2 Om(1,2)^2 Om(1,3)^2; Om(2,1)^2 Om(2,2)^2 Om(2,3)^2; Om(3,1)^2 Om(3,2)^2 Om(3,3)^2];
K2 = [Om(1,2)*Om(1,3) Om(1,3)*Om(1,1) Om(1,1)*Om(1,2); Om(2,2)*Om(2,3) Om(2,3)*Om(2,1) Om(2,1)*Om(2,2); Om(3,2)*Om(3,3) Om(3,3)*Om(3,1) Om(3,1)*Om(3,2)];
K3 = [Om(2,1)*Om(3,1) Om(2,2)*Om(3,2) Om(2,3)*Om(3,3); Om(3,1)*Om(1,1) Om(3,2)*Om(1,2) Om(3,3)*Om(1,3); Om(1,1)*Om(2,1) Om(1,2)*Om(2,2) Om(1,3)*Om(2,3)];
K4 = [Om(2,2)*Om(3,3)+Om(2,3)*Om(3,2) Om(2,3)*Om(3,1)+Om(2,1)*Om(3,3) Om(2,1)*Om(3,2)+Om(2,2)*Om(3,1); Om(3,2)*Om(1,3)+Om(3,3)*Om(1,2) Om(3,3)*Om(1,1)+Om(3,1)*Om(1,3) Om(3,1)*Om(1,2)+Om(3,2)*Om(1,1); Om(1,2)*Om(2,3)+Om(1,3)*Om(2,2) Om(1,3)*Om(2,1)+Om(1,1)*Om(2,3) Om(1,1)*Om(2,2)+Om(1,2)*Om(2,1)];
K = [K1 2*K2;K3 K4];
% K = [m*m n*n 0 0 0 2*m*n; n*n m*m 0 0 0 -2*m*n; 0 0 1 0 0 0; 0 0 0 m -n 0; 0 0 0 n m 0; -m*n m*n 0 0 0 m*m-n*n];
Kt=K.';
Cv = K*C*Kt;
Cv_GPa = Cv*1e-3

S=inv(Cv)
disp('Ex Ey Ez Gyz Gzx Gxy vzx vxz');
eng = [1e-3/S(1,1) 1e-3/S(2,2) 1e-3/S(3,3) 1e-3/S(4,4) 1e-3/S(5,5) 1e-3/S(6,6) -S(1,3)/S(3,3) -S(3,1)/S(1,1)];
disp(eng);

% Find average isotropic material Properties
% Voigt Average - uniform strain, over-estimate of stresses (Get from C)
% bulk modulus
Bv=0;
for i=1:3
Bv=Bv+sum(Cv(i,1:3));
end
Bv=Bv/9;
% shear modulus
Gv=(Cv(1,1)+Cv(2,2)+Cv(3,3))/15-(Cv(1,2)+Cv(2,3)+Cv(1,3))/15+(Cv(4,4)+Cv(5,5)+Cv(6,6))/5;
% Youngs Modulus & Poissons Ratio
Ev=(9*Bv*Gv)/(3*Bv+Gv);
nuv=(3*Bv-2*Gv)/2/(3*Bv+Gv);
% Reuss Average - uniform stres, under-estimate of stressses (Get from S)
Br=0;
for i=1:3
Br=Br+sum(S(i,1:3));
end
Br=1/Br;
% Shear Modulus
Gr=4/15*(S(1,1)+S(2,2)+S(3,3))-4/15*(S(1,2)+S(2,3)+S(1,3))+1/5*(S(4,4)+S(5,5)+S(6,6));
Gr=1/Gr;
% Youngs Modulus & Poissons Ratio
Er=(9*Br*Gr)/(3*Br+Gr);
nur=(3*Br-2*Gr)/2/(3*Br+Gr);
fprintf('\nISOTROPIC PROPERTIES:\n')
fprintf(' Voigt Avg Properties: Ev=%6.2f, nuv=%6.2f, Gv=%6.2f , Bv=%6.2f,  kv(G/B)=%6.2f \n',Ev,nuv,Gv,Bv,Gv/Bv)
fprintf(' Reuss Avg Properties: Er=%6.2f, nur=%6.2f, Gr=%6.2f , Br=%6.2f,  kr(G/B)=%6.2f \n',Er,nur,Gr,Br,Gr/Br)
%--------------------------------------------------------------------------
% get fourth order C
C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                m=0;
                n=0;
                if i==j
                    m=i;
                elseif ((i==2 && j==3)||(i==3 && j==2))
                    m=4;
                elseif ((i==1 && j==3)||(i==3 && j==1))
                    m=5;
                elseif ((i==1 && j==2)||(i==2 && j==1))
                    m=6;
                end
                if k==l
                    n=l;
                elseif ((k==2 && l==3)||(k==3 && l==2))
                    n=4;
                elseif ((k==1 && l==3)||(k==3 && l==1))
                    n=5;
                elseif ((k==1 && l==2)||(k==2 && l==1))
                    n=6;
                end
                C(i,j,k,l)=Cv(m,n);
            end
        end
    end
end
%----CALCULATE KK----------------------------------------------------------
%permutation or Levi-Cevita 3x3x3 tensor (using linear indexing of matrix)
lc = zeros(3,3,3);
lc([8 12 22]) = 1;
lc([6 16 20]) = -1;
% Integral for equation 2.2
t=line_dir/norm(line_dir); %disocation line direction -- eq 2.3
r=norm(t);
theta=atan2(t(2),t(1));
phi=acos(t(3)/r);
a=[sin(theta) -cos(theta) 0]; %eq 2.5
d=[cos(phi)*cos(theta) cos(phi)*sin(theta) -sin(phi)]; % eq 2.5
%Integral Range: psi ranges 0 to pi
psi_1=0;
psi_2=pi;
Ninc=200; % # of midpoint summation intervals
psi_inc=(psi_2-psi_1)/Ninc;
intzM=zeros(3,3,3,3);
intzMiso=zeros(3,3,3,3);
LAM3=zeros(3,3);
for inc=1:Ninc %midpoint integration loop for intgral of christophel matrix
    psi=psi_1+(inc-1)*psi_inc;
    z(1:3) =a*cos(psi)+d*sin(psi); %eq 2.4
    dz(1:3)=-a*sin(psi)+d*cos(psi); %eq 2.4
    % Christoffel stiffness -- C_ijrs z_r z_s
    M=zeros(3,3);
    for i=1:3
        for r=1:3
            for j=1:3
                for s=1:3
                    M(i,r)=M(i,r)+C(i,j,r,s)*z(j)*z(s);
                end
            end
        end
    end
    %midpoint integration of Christophel terms for eq 2.2
    invM=inv(M);
    for s=1:3
        for n=1:3
            for i=1:3
                for r=1:3
                    intzM(s,n,i,r)=intzM(s,n,i,r) + (z(s)*dz(n)*invM(i,r))*psi_inc;
                end
            end
        end
    end
end
% SOLVE FOR ANISOTROPIC DISLOCATION COMPLIANCE OR ORIENTATION PREFACTOR TENSORS
% this is the inverse of LAMBDA given by Sun & Beltz (JMPS 1994) etc.
KK_BA=zeros(3,3); % Barnett & Asaro JMPS v.20 1972
for m=1:3
    for g=1:3
        for p=1:3
            for j=1:3
                for w=1:3
                    for n=1:3
                        for i=1:3
                            for r=1:3
                                for s=1:3
                                    % -- Equation 2.2
                                    KK_BA(m,g)=KK_BA(m,g)+1/8/pi^2*lc(p,j,w)*t(j)...
                                        *(C(n,g,i,p)*C(w,m,r,s)+C(n,m,i,p)*C(w,g,r,s))...
                                        *intzM(s,n,i,r);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%%%%---OUTPUT RESULTS---
%
% Dislocation energy
Eaniso=bvec'*(KK_BA)*bvec * 1e6 * 1e9; % eq 2.1
cos_theta=dot(bvec,line_dir/norm(line_dir))/norm(bvec);
Eiso_v=1/4/pi*Gv*(1-nuv*cos_theta^2)/(1-nuv)*(norm(bvec))^2*1e6*1e9;%eq 2.6
Eiso_r=1/4/pi*Gr*(1-nur*cos_theta^2)/(1-nur)*(norm(bvec))^2*1e6*1e9;

chosen
fprintf('\nENERGY PRELOG FACTOR (Jm^-1*1e-9): \n\n')
fprintf(' Anisotropic = %-6g\n',Eaniso)
fprintf(' Isotropic Voigt = %-6g\n',Eiso_v)
fprintf(' Isotropic Reuss = %-6g\n',Eiso_r)
fprintf(' b^2 = %-6g\n',1e20*norm(bvec)^2)
fprintf('%-6g %-6g %-6g\n',Eaniso,1e20*norm(bvec)^2,Eiso_v)
%
% splitting distance between two edge partials
%d_aniso=(4*pi*bvec'*(KK_BA)*bvec) /2/pi/(y_sf*1e-6)*1e10; %eq 3.17
%d_isov= (Gv/(1-nuv))/2/pi/(y_sf*1e-6)*(norm(bvec))^2*1e10;
%d_isor= (Gv/(1-nur))/2/pi/(y_sf*1e-6)*(norm(bvec))^2*1e10;
%fprintf('\n\nEDGE PARTIAL SPLITTING DISTANCE (A): \n\n')
%fprintf(' Anisotropic = %-6g\n',d_aniso)
%fprintf(' Isotropic Voigt = %-6g\n',d_isov)
%fprintf(' Isotropic Reuss = %-6g\n',d_isor)
