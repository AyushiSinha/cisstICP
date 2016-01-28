clear, clc

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

% load samples
[X] = PointsLoad(samplesPath);

% load initial offset
fid = fopen(offsetXfmPath);
if fid == -1
  error(['failed to open file: ',offsetXfmPath]);
end
c = textscan(fid,'%f',12);
fclose(fid);
c = c{:};
Ri = reshape(c(1:9),3,3)';
% normalize rotation
[u, s, v] = svd(Ri);
Ri = u*v';
ti = c(10:12);
Fi = getFrm3(Ri,ti);

disp('Initial Offset')
disp(Fi)

% set ICP options
opt = objOptICP();
opt.term_dAng = 0.01*pi/180;
opt.term_dPos = 0.01;
opt.term_holdIter = 2;
opt.maxIter = 200;

% set ICP algorithm
flags.bEnableDebug = false;
flags.bTargetAsMesh = true;
alg = algICP_StdICP();
alg.Initialize( V,T,Tn, X, flags, {mexDir,mexInterfaceDir} );

%--- Run ICP ---%
extrasICP.bEnableDebug = false;

[Freg] = IterateICP( alg, opt, Fi, extrasICP );
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
