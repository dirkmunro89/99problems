%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%function top88(nelx,nely,volfrac,penal,rmin,ft)
mmatopinit;
nelx=60;
nely=60;
volfrac=0.1;
penal=3.0;
rmin=3.3;
ft=2;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1*0.,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)-1]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(0.1,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
fscl=1;
%% START ITERATION
while change > 0.001 && loop < 2000
  loop = loop + 1;
  %% SELF-WEIGHT LOAD
  Fx=F;
  gv = -1.;
  bf = 1./4.*[0. gv 0. gv 0. gv 0. gv]';
  for e = 1:nelx*nely
    Fx(edofMat(e,:)) = Fx(edofMat(e,:)) + bf*xPena(0.,1.,xPhys(e));
  end
  %% FE-ANALYSIS
  xPhyp = xPena(0.0,penal,xPhys);
  dxPhyp = dxPena(0.0,penal,xPhys);
  sK = reshape(KE(:)*(Emin+xPhyp(:)'*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\Fx(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  f = sum(sum((Emin+xPhyp*(E0-Emin)).*ce));
  df = -(E0-Emin)*dxPhyp.*ce;
  for e = 1:nelx*nely
    df(e) = df(e) + 2.*bf'*U(edofMat(e,:))*dxPena(0.,1.,xPhys(e));
  end
  dv = -ones(nely,nelx)/volfrac/nelx/nely;
  %% SCALING
  if loop == 1
    fscl = f;
  end
  f=f/fscl;
  df=df/fscl;
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    df(:) = H*(x(:).*df(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    df(:) = H*(df(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% DATA STRUCTURES FOR THE MMA
  f0val=f;
  df0dx=reshape(df,[n,1]);
  fval=[-sum(xPhys(:))/volfrac/nelx/nely + 1];
  dfdx = [reshape(dv,[1,n])];
  %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
  %% HISTORY AND CHANGE TO TOPOPT STRUCTURE
  xold2 = xold1;
  xold1 = xval;
  xval = xmma;
  xnew = reshape(xval,[nely,nelx]);
  %% FILTERING
  if ft == 1
    xPhys = xnew;
  elseif ft == 2
    xPhys(:) = (H*xnew(:))./Hs;
  end
  %% CHANGE AND OUTPUT
  change = max(abs(xnew(:)-x(:)));
  x=xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4e Vol.:%7.3e ch.:%7.3f\n',loop,f*fscl, ...
    fval,change);
  %% PLOT DENSITIES
% fig = figure('visible', 'off');
% colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off;
% saveas(fig,'topo','png');
% close(fig);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

