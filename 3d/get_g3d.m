function [g,Dg,H,actual_vols] = get_g3d(w,X,target_vols,bx,periodic,varargin)

% [g,Dg,H,actual_vols] = get_g(w,X,target_vols,bx,periodic,OPTcheckneighbours)
%
% WARNING! By selecting periodic=true and checkneighbours=false there is a risk that the Hessian is incorrect
% Proceed at your own risk! If all cell volumes are equal then for N greater than around 50 in a cube the result is correct almost all the time
%
% Input arguments
%         w            - is the weights, an Nx1 array
%         X            - is the positions, an Nx3 array
%         target_vols  - target volumes, an Nx1 array
%         bx           - box size, a 1x3 array
%         periodic     - periodic flag (a boolean true/false to indicate periodic or not)
%         OPTcheckneighbours - optional argument to check neighbours in the
%         periodic case, default is TRUE. This argument is ignored if
%         periodic=false
%    
% Return arguments
%         g            - function g(w;x)
%         Dg           - gradient of g wrt w
%         H            - the Hessian d^2g, a sparse matrix
%         actual_vols  - the actual volumes of the Laguerre diagram with seeds X and weights w
    
    %% Catch errors and process the arguments
    
    if(nargin==6)
        checkneighbours=varargin{1};
    else
        checkneighbours=true;
    end
    
    % If we are periodic then we call get_g_voro, by default checkneighbours is true
    if(periodic)
        
        % We must map the seeds into the primary domain
        X=mod(X*diag(1./bx),1)*diag(bx);
        
        [g,Dg,H,actual_vols,success]=get_g_voro(w,X,target_vols,bx,periodic,checkneighbours);
        % If get_g_voro detected a problem then it exits with success=false
        % In this case we then call get_g_kper, as this is robust (but slow)
        if(~success)
            warning('Neighbour problem detected, switching to get_g_kper, may be slower');
            [g,Dg,H,actual_vols]=get_g_kper(w,X,target_vols,bx);
        end
    else
        % If we are non-periodic then get_g_voro is robust
        [g,Dg,H,actual_vols,~]=get_g_voro(w,X,target_vols,bx,periodic);
    end
end

