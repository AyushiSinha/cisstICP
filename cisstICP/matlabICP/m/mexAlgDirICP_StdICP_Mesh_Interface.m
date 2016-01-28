classdef mexAlgDirICP_StdICP_Mesh_Interface < handle
% MATLAB class wrapper to an underlying C++ class

  properties (SetAccess = private)
    
    h;  % Handle to the underlying C++ class instance
    
  end
  
  methods
    
    %% Constructor - Create a new C++ class instance 
    function this = mexAlgDirICP_StdICP_Mesh_Interface()
      this.h = mexAlgDirICP_StdICP_Mesh('new');
    end

    %% Destructor - Destroy the C++ class instance
    function delete(this)
      mexAlgDirICP_StdICP_Mesh('delete', this.h);
    end
    
    %% Initialize Algorithm
    % inputs:
    %   V     ~ target vertices   Nv x 3
    %   T,Tn  ~ target triangles and triangle normals   Nt x 3
    %   Xp    ~ source points             Ns x 3
    %   Xn    ~ source normals            Ns x 3
    function Initialize(this, V,T,Tn, Xp,Xn)
      % convert triangles to C++ indexing (base-0) and to integer data type
      Tcpp = int32(T - ones(size(T)));
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  mesh
      % //    V ~ 3D vertex positions     (3 x Nv double)
      % //    T ~ triangle vertex indices (3 x Nt integer)
      % //    N ~ 3D triangle normals     (3 x Nt double)
      % //  sample points                 (3 x Ns double)
      % //  sample normals                (3 x Ns double)
      % //
      mexAlgDirICP_StdICP_Mesh('Initialize',this.h,...
        V',Tcpp',Tn',...
        Xp',Xn');
    end
    
    %% Set Samples
    % inputs:
    %   Xp    ~ source points             Ns x 3
    %   Xn    ~ source normals            Ns x 3
    function SetSamples(this, Xp,Xn )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  sample points                 (3 x Ns double)
      % //  sample normals                (3 x Ns double)
      % //
      mexAlgDirICP_StdICP_Mesh('SetSamples',this.h,...
        Xp',Xn' );
    end
        
    %% ICP: Initialize Parameters
    % inputs:
    %   F0   ~ initial registration guess   (4 x 4 double)
    function ICP_InitializeParameters(this, F0)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % //  Freg ~ [R,t]      (3 x 4 double)
      mexAlgDirICP_StdICP_Mesh('ICP_InitializeParameters',this.h, F0(1:3,:));
    end    
        
    %% ICP: Compute Matches
    % computes matches and post-match parameters
    %
    % outputs:
    %   Yp   ~ match positions    N x 3
    %   Yn   ~ match orientations N x 3
    function [Yp, Yn] = ICP_ComputeMatches(this)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  SampleMatchPts ~ matches for 3D sample pts    (3 x numPts double)  [optional]      
      [Yp, Yn] = mexAlgDirICP_StdICP_Mesh('ICP_ComputeMatches', this.h);
      Yp = Yp';
      Yn = Yn';
    end    
    
    %% ICP: Register Matches
    % registers matches and computes post-registration parameters
    %
    % outputs:
    %   F   ~ registration parameters
    function [F] = ICP_RegisterMatches(this)
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  Freg ~ [R,t]      (3 x 4 double)    (Y = Freg * X)
      % //
      Freg = mexAlgDirICP_StdICP_Mesh('ICP_RegisterMatches', this.h);
      F = getFrm3(Freg(1:3,1:3),Freg(1:3,4));
    end        

    function fval = ICP_EvaluateErrorFunction( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  fval  ~ value of error function (double)
      % //      
      fval = mexAlgDirICP_StdICP_Mesh('ICP_EvaluateErrorFunction', this.h);
    end
    
    function bTerminate = ICP_Terminate( this )
      % // Expected Input:
      % //  cmd
      % //  class handle
      % // 
      % // Output:
      % //  bTerminate  ~ boolean specifying algorithm-specific termination (logical scalar)
      % //
      bTerminate = mexAlgDirICP_StdICP_Mesh('ICP_Terminate', this.h);
    end
  end
end