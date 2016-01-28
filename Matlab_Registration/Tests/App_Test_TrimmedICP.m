clear, clc

scale0 = 1.0;
bEstimateScale = true;
addOutliersRate = 0.1;
trimRate = 0.1;
nSamps = 100;

inputDir = './TestData_StdICP';

meshPath = [inputDir,'/','SaveMesh.mesh'];
samplesPath = [inputDir,'/','SaveNoisySamples.pts'];
offsetXfmPath = [inputDir,'/','SaveOffsetXfm.txt'];

%mexDir = 'D:\Code\Repos_Git\RegistrationTools\matlabICP_build_x64\Release';
mexDir = 'D:\Code\Repos_Git\RegistrationTools\matlabICP_build_x64_debug\Debug';
mexInterfaceDir = 'D:\Code\Repos_Git\RegistrationTools\cisstICP\matlabICP\m';
%mexPath = [mexDir,'/','mexAlgICP_StdICP_Mesh.mexw64'];

% load mesh
[V,T,Tn] = MeshLoad(meshPath);

% generate samples
sampsSeed = 11;
[X] = GenerateRandomMeshSamples(V,T,Tn, nSamps, sampsSeed);

% load samples
%[X] = PointsLoad(samplesPath);

% generate outlier samples
nOutliers = size(X,1)*addOutliersRate;
outlierSeed = 5;
outlierDist = 20;
[oXp, oXn] = GenerateRandomMeshSamples(V,T,Tn, nOutliers, outlierSeed);
oXp = oXp + oXn * 20;
X = [X;oXp];

% generate initial offset
offsetT = [5,5.1];
offsetR = [5,5.1];
offsetSeed = 7;
Fi = GenerateRandomTransforms(offsetT(1),offsetT(2),offsetR(1),offsetR(2),1,offsetSeed);

% % load initial offset
% fid = fopen(offsetXfmPath);
% if fid == -1
%   error(['failed to open file: ',offsetXfmPath]);
% end
% c = textscan(fid,'%f',12);
% fclose(fid);
% c = c{:};
% Ri = reshape(c(1:9),3,3)';
% ti = c(10:12);
% Fi = getFrm3(Ri,ti,true);

disp('Initial Offset')
disp(Fi)

% set ICP options
opt = objOptICP();
opt.term_dAng = 0.01*pi/180;
opt.term_dPos = 0.01;
opt.term_holdIter = 2;
opt.maxIter = 100;

% set ICP algorithm
flags.bEnableDebug = false;
flags.bTargetAsMesh = true;
flags.bEstimateScale = bEstimateScale;
flags.bPlotIterations = true;
flags.trimRate = trimRate;
algPayload.scale0 = scale0;

alg = algICP_TrimmedICP();
alg.Initialize( V,T,Tn, X, flags, {mexDir,mexInterfaceDir} );

%--- Run ICP ---%
extrasICP.bEnableDebug = false;
extrasICP.bPlotIter = false;
extrasICP.Vplot = V;
extrasICP.Tplot = T;
extrasICP.Xplot = X;
extrasICP.XplotGT = [];

[Freg] = IterateICP( alg, opt, Fi, extrasICP, algPayload );
disp('Freg')
disp(Freg)

% compute registraiton error
% assuming ground truth is identity
[~,AngErr] = rot2AxisAngle(getRot(Freg));
PosErr = norm(getPos(Freg));

% compare to initial offset
[~,AngErrInit] = rot2AxisAngle(getRot(Fi));
PosErrInit = norm(getPos(Fi));

disp('AngErrInit  PosErrInit')
disp([AngErrInit*180/pi PosErrInit])

disp('AngErr      PosErr')
disp([AngErr*180/pi PosErr])

% cleanup
alg.delete();   % calling "clear" also does the trick
