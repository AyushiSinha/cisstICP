clear, clc

nSamps = 100;

seedSamps = 7;
seedXfms = 11;
seedL = 3;

offsetTrans = [3, 3.01];
offsetRot = [3, 3.01];

inputDir = './TestData_StdICP';

meshPath = [inputDir,'/','SaveMesh.mesh'];
%samplesPath = [inputDir,'/','SaveNoisySamples.pts'];
%offsetXfmPath = [inputDir,'/','SaveOffsetXfm.txt'];

%mexDir = 'D:\Code\Repos_Git\RegistrationTools\matlabICP_build_x64\Release';
mexDir = 'D:\Code\Repos_Git\RegistrationTools\matlabICP_build_x64_debug\Debug';
mexInterfaceDir = 'D:\Code\Repos_Git\RegistrationTools\cisstICP\matlabICP\m';
%mexPath = [mexDir,'/','mexAlgICP_StdICP_Mesh.mexw64'];

% load mesh
[V,T,Tn] = MeshLoad(meshPath);

% generate random samples
[XpGt, XnGt, Xdatum] = GenerateRandomMeshSamples(V,T,Tn,nSamps,seedSamps);

% generate misalignment
[Foffset] = GenerateRandomTransforms(...
  offsetTrans(1),offsetTrans(2),offsetRot(1),offsetRot(2),...
  1, seedXfms);
Fgt = invFrm3(Foffset);

% apply misalignment
Xp = bsxfun(@plus, XpGt * getRot(Foffset)', getPos(Foffset)');
Xn = XnGt * getRot(Foffset)';

% define noise model
M = repmat(eye(3,3),1,1,nSamps);
k = ones(nSamps,1)*100;
B = zeros(nSamps,1);
L = GenerateRandom_KentMajorMinorAxis(Xn,seedL);

% define initial guess
Ri = eye(3);
ti = zeros(3,1);
Fi = getFrm3(Ri,ti);

% set ICP options
opt = objOptSettingsICP();
opt.term_dAng = 0.01*pi/180;
opt.term_dPos = 0.01;
opt.term_holdIter = 2;
opt.maxIter = 100;



%--- Kent ---%

% set ICP algorithm
bEnableDebug = 0;
bTargetAsMesh = 1;
alg = algDirICP_Kent();
alg.Initialize(...
  V,T,Tn,...
  Xp,Xn,M,k,B,L,...
  bTargetAsMesh,bEnableDebug, {mexDir,mexInterfaceDir} );

%--- Run ICP ---%
[Freg] = IterateICP( alg, opt, Fi, bEnableDebug );
disp('Freg')
disp(Freg)

% compute registration frame error
dF = Fgt*invFrm3(Freg);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));
% % assuming ground truth is identity
% [~,AngErr] = rot2AxisAngle(getRot(Freg));
% PosErr = norm(getPos(Freg));

% compare to initial offset
[~,AngErrInit] = rot2AxisAngle(getRot(Foffset));
PosErrInit = norm(getPos(Foffset));

disp('AngErrInit  PosErrInit')
disp([AngErrInit*180/pi PosErrInit])

disp('AngErr      PosErr')
disp([AngErr*180/pi PosErr])

% cleanup
alg.delete();   % calling "clear" also does the trick



%--- Standard ICP ---%

% set ICP algorithm
bEnableDebug = 0;
bTargetAsMesh = 1;
alg = algDirICP_StdICP();
alg.Initialize(...
  V,T,Tn,...
  Xp,Xn,...
  bTargetAsMesh,bEnableDebug, {mexDir,mexInterfaceDir} );

%--- Run ICP ---%
[Freg2] = IterateICP( alg, opt, Fi, bEnableDebug );
disp('Freg')
disp(Freg2)

% compute registration frame error
dF = Fgt*invFrm3(Freg2);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));

% compare to initial offset
[~,AngErrInit] = rot2AxisAngle(getRot(Foffset));
PosErrInit = norm(getPos(Foffset));

disp('AngErrInit  PosErrInit')
disp([AngErrInit*180/pi PosErrInit])

disp('AngErr      PosErr')
disp([AngErr*180/pi PosErr])

% cleanup
alg.delete();   % calling "clear" also does the trick



%--- Fisher ---%

k_init = k(1);
sigma2_init = 1.0;
wRpos = 0.5;
bDynamicEst = false;

% set ICP algorithm
bEnableDebug = 0;
bTargetAsMesh = 1;
alg = algDirICP_vMFG();
alg.Initialize(...
  V,T,Tn,...
  Xp,Xn, k_init,sigma2_init,wRpos,bDynamicEst,...
  ...%Xp,Xn, [],[],[],[],...  
  bTargetAsMesh,bEnableDebug, {mexDir,mexInterfaceDir} );

%--- Run ICP ---%
[Freg3] = IterateICP( alg, opt, Fi, bEnableDebug );
disp('Freg')
disp(Freg3)

% compute registration frame error
dF = Fgt*invFrm3(Freg3);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));

% compare to initial offset
[~,AngErrInit] = rot2AxisAngle(getRot(Foffset));
PosErrInit = norm(getPos(Foffset));

disp('AngErrInit  PosErrInit')
disp([AngErrInit*180/pi PosErrInit])

disp('AngErr      PosErr')
disp([AngErr*180/pi PosErr])

% cleanup
alg.delete();   % calling "clear" also does the trick


% =============================

disp(' ')
disp(' ');

% Compare Kent and StdICP
dF = Freg * invFrm3(Freg2);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));

disp('Kent vs StdICP:')
disp('dAng      dPos')
disp([AngErr*180/pi PosErr])


% Compare Kent and Fisher
dF = Freg * invFrm3(Freg3);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));

disp('Kent vs Fisher:')
disp('dAng      dPos')
disp([AngErr*180/pi PosErr])


% load samples
%[X] = PointsLoad(samplesPath);
% % load initial offset
% fid = fopen(offsetXfmPath);
% if fid == -1
%   error(['failed to open file: ',offsetXfmPath]);
% end
% c = textscan(fid,'%f',12);
% fclose(fid);
% c = c{:};
% Ri = reshape(c(1:9),3,3)';
% % normalize rotation
% [u, s, v] = svd(Ri);
% Ri = u*v';
% ti = c(10:12);
% Fi = getFrm3(Ri,ti);
% 
% disp('Initial Offset')
% disp(Fi)