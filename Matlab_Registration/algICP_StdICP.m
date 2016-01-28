classdef algICP_StdICP < algICP

properties

  bEnableDebug;
  bTargetAsMesh;
  
  % sample data
  X;   % source pts        N x 3
  Y;   % matches 3D pts    N x 3
  
  % target shape
  % mesh and point cloud target:
  V;   % vertices            Nv x 3
  % mesh target only:
  T;   % triangle faces (vertex indices)         Nt x 3
  Tn;  % triangle normals                        Nt x 3
  
  % outliers
  outliers; % logical array of outliers
  
  % member objects
  mexAlgStdICP;   % mex object for StdICP algorithm
end

methods
  
  %% Constructor
  function [this] = algICP_StdICP()
  end
  
  %% Delete
  function delete( this )
    % delete mex objects
    this.mexAlgStdICP.delete();
  end  
  
  %% Initialize Algorithm
  % V       ~ target vertices     Nv x 3
  % T,Tn    ~ target triangles and triangle normals (for mesh target only)  Nt x 3
  % X       ~ source points       N x 3  
  % flags
  %   bEnableDebug    ~ enable debug messages
  %   bTargetAsMesh   ~ set point cloud or mesh target
  % paths   ~ include paths (i.e. mex file path)  (cell array)
  function Initialize( this, ...
    V,T,Tn, X, flags, paths )
    
    % initialize superclass
    Initialize@algICP( this );
        
    % debug enable
    if ~exist('flags','var') || ~isfield(flags,'bEnableDebug') || isempty(flags.bEnableDebug)
      flags.bEnableDebug = 0;
    end
    this.bEnableDebug = flags.bEnableDebug;
    
    % source shape
    if ~exist('X','var') || isempty(X)
      %error('Missing source shape input X')
      warning('Missing source shape input X')
      this.X = [];
    else      
      if size(X,2) ~= 3
        error('Invalid size of source shape (X)')
      end
      this.X = X;
    end
        
    % target shape
    if ~exist('flags','var') || ~isfield(flags,'bTargetAsMesh') || isempty(flags.bTargetAsMesh)
      flags.bTargetAsMesh = 1;
    end
    if size(V,2) ~= 3
      error('Invalid size of target vertices (V)')
    end
    this.bTargetAsMesh = flags.bTargetAsMesh;
    this.V = V;
    if (this.bTargetAsMesh)
      this.T = T;
      this.Tn = Tn;
    else
      this.T = [];
      this.Tn = [];
    end
    
    % include paths
    %  (i.e. for directories containing mex libraries)
    for i=1:length(paths)
      addpath( paths{i} );
    end
        
    % Mex objects
    if (this.bTargetAsMesh)
      utlDebugMsg( this.bEnableDebug,'Creating mex object: mexAlgICP_StdICP_Mesh...\n')
      this.mexAlgStdICP = mexAlgICP_StdICP_Mesh_Interface();
      utlDebugMsg( this.bEnableDebug,' ...done\n');
      
      utlDebugMsg( this.bEnableDebug,'Initializing mexAlgStdICP...\n')
      this.mexAlgStdICP.Initialize( this.X, this.V,this.T,this.Tn );
      utlDebugMsg( this.bEnableDebug,' ...done\n');
    else
      error('mex not supported: mexAlgStdICP_PointCloud');
    end    
  end  
  
  %% Set Samples
  % X  ~ source points       N x 3
  function SetSamples( this, X )  
    this.X = X;
    this.mexAlgStdICP.SetSamples( X );
  end
  
  %% ICP: Initialize Registration
  function ICP_InitializeParameters( this, F0, algPayload )       
    
    % initial guess
    if ~exist('F0','var') || isempty(F0)
      F0 = getFrm3(eye(3),zeros(3,1));
    end
    
    % initialize superclass
    ICP_InitializeParameters@algICP(this, F0);
    
    utlDebugMsg( this.bEnableDebug,'initializing ICP parameters...\n')
    this.mexAlgStdICP.ICP_InitializeParameters( F0 );
    utlDebugMsg( this.bEnableDebug,' ...done\n');
  end 
  
  %% ICP: Compute Matches
  function [ Y ] = ICP_ComputeMatches( this )
    % compute matches
    %  also computes post-match parameters
    utlDebugMsg( this.bEnableDebug,'computing matches...\n' );
    this.Y = this.mexAlgStdICP.ICP_ComputeMatches();
    Y = this.Y;
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
  end

  %% ICP: Register Matches
  function [F,dRrod,dt] = ICP_RegisterMatches( this )
   
    % save prior registration
    Fprior = this.F;
    
    % register matches
    %  also computes post-registration parameters
    utlDebugMsg( this.bEnableDebug,'registering matches...\n' );
    [this.F] = this.mexAlgStdICP.ICP_RegisterMatches();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );    
    
    % compute change in registration
    dF = this.F * invFrm3(Fprior);
    this.dRrod = rot2rodrigues(getRot(dF));
    this.dt = getPos(dF);
    
    % return values
    F = this.F;
    dRrod = this.dRrod;
    dt = this.dt;
  end
  
  %% ICP: Update Parameters (Post Register)
  %   use this function when computing parameters following registration
  function ICP_UpdateParameters_PostRegister( this, Freg )
    % update post-registration parameters
    utlDebugMsg( this.bEnableDebug,'updating post-registration parameters...\n' );
    this.mexAlgStdICP.ICP_UpdateParameters_PostRegister( Freg );
    utlDebugMsg( this.bEnableDebug,' ...done\n' );                
  end
  
  %% ICP: Evaluate Error Function
  function fval = ICP_EvaluateErrorFunction( this )
    utlDebugMsg( this.bEnableDebug,'evaluating error function...\n' );
    this.fval = this.mexAlgStdICP.ICP_EvaluateErrorFunction();
    utlDebugMsg( this.bEnableDebug,' ...done\n' );
    
    % return values
    fval = this.fval;
  end
  
end

end
