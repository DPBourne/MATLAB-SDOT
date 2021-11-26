function [g,Dg,H,actual_vols,success] = get_g_voro(w,X,target_vols,bx,periodic,varargin)
        
    %% WARNING! Is only accurate when the boundary between two cells is composed of a single facet
    %
    % [g,Dg,H,actual_vols,success] = get_g(w,X,target_vols,bx,periodic,OPTcheckneighbours)
    %
    % Input arguments
    %         w            - is the weights, an Nx1 array
    %         X            - is the positions, an Nx3 array
    %         target_vols  - target volumes, an Nx1 array
    %         bx           - box size, a 1x3 array
    %         periodic     - periodic flag (a boolean true/false to indicate periodic or not)
    %         OPTcheckneighbours - is an optional boolean true/false to indicate whether the function checks for problems with neighbours (multiple facets or self facets)
    
    % Return arguments
    %         g            - function g(w;x)
    %         Dg           - gradient of g wrt w
    %         H            - the Hessian d^2g, a sparse matrix
    %         actual_vols  - the actual volumes of the Laguerre diagram with seeds X and weights w
    %         success      - a flag to indicate everything has been calculated successfully
    
    %% Catch errors and process arguments
    
    % Here we check the input arguments to see if the optional 'checkneighbours' flag has been set
    if(nargin==6)
        checkneighbours=varargin{1};
    else
        checkneighbours=false;
    end
    
    [Nw,Mw]=size(w);    
    
    if(Mw~=1)
        error('w should be an N x 1 array where N is the number of cells');
    end

    [NX,MX]=size(X);

    if(MX~=3)
        error('X should be an N x 3 array where N is the number of cells');
    end

    if(NX~=Nw)
        error('The number of cells represented by w and X disagree, w: %d and X %d',Nw,NX);
    end

    [Ntv,Mtv]=size(target_vols);
    if(Mtv~=1)
        error('target_vols should be an N x 1 array where N is the number of cells');
    end

    if(Ntv~=Nw)
        error('The number of cells represented by w and target_vols disagree, w: %d and target_vols: %d',Nw,Ntv);
    end

    total_vol=sum(target_vols);
    bx_vol=bx(1)*bx(2)*bx(3);

    if(abs(total_vol-bx_vol)/bx_vol>1e-10)
        error('The sum of the target volumes is different to the volume of the box');
    end

    
    %% Computations
    
    % Use mexPDallfaces to calculate the Laguerre diagram
    [actual_vols,transport_costs,~,vfn]=mexPDallfaces(bx,X,w,periodic);
    
    % Gradient of g
    Dg=actual_vols-target_vols;

    % Definition of g, a convex function of the weights 
    g=dot(Dg,w)-sum(transport_costs); 
    
    % Preallocate memory for possible speed improvement
    % Upper bound on number of off-diagonal nonzero entries in the Hessian is the sum of the number of neighbours over all cells 
    numNZ=sum(cellfun(@length,vfn(:,3)));
    % Make a zero array of the correct size
    Z=zeros(1,numNZ);
    
    % Indices of matrix to store inter-seed distances and face areas
    AD_spvals_i=Z;
    AD_spvals_j=Z;
    
    % Values
    D_spvals_val=Z;
    A_spvals_val=Z;
    
    % Index to keep track of the number of entries in the sparse matrix
    r=1;
    
    % Loop over cells
    
    for i=1:Nw
        % Coordinates of the seed location of the ith cell
        xi=X(i,:);

        % The vertices for the ith cell
        Vertices_Cell_i=vfn{i,1};
        
        % The list of neighbour indices of the ith cell - boundaries have negative indices
        N_i=vfn{i,3};
        
        % The vertex indices for each of the faces between the cell and the neighbours in N_i
        Face_Index_i=vfn{i,2}(:,1);
        Face_Areas_i=vfn{i,2}(:,2);
        Face_Normals_i=vfn{i,2}(:,3);
        
        % The number of neighbours
        Num_N_i = size(N_i);
        
        % ERROR CHECKING, to see if a cell has a boundary with another cell composed of multiple facets, or a boundary between the cell and a periodic copy of itself
        
        if(checkneighbours)
            if(length(N_i)~=length(unique(N_i)))
                %error(sprintf('Cell with index %d has a boundary with another cell composed of multiple facets',i));
                success=false;
                H=[];
                return;
            end
            
            if(ismember(i,N_i))
                %error(sprintf('Cell with index %d has part of its boundary arising from an intersection with a periodic copy of itself',i));
                success=false;
                H=[];
                return;
            end
        end
        
        % Calculate H_ij for j in N_i 
        for j=1:Num_N_i
            
            % The cell index of neighbour j
            k = N_i(j);
            
            % We only need to calculate if the neighbour is a not a part of the boundary (non-positive k)
            if(k>0)
                
                xk=X(k,:);
                
                % The indices of the vertices for the face between cell i and cell k (neighbour j)
                Face_Vertex_Indices=Face_Index_i{j};
                
                % The face vertices for the face between cell i and cell k (neighbour j)
                Face_Vertices=Vertices_Cell_i(Face_Vertex_Indices,:);
                
                % Area of the face 
                a_ik=Face_Areas_i{j};
                
                AD_spvals_i(r)=i;
                AD_spvals_j(r)=k;
                
                A_spvals_val(r)=a_ik;
                
                % %%%%%%%%%%%%%%%%% For checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                a_ik_polyarea = polygonArea3d(Face_Vertices);
                %                if(abs(a_ik-a_ik_polyarea)>1e-10)
                %                    disp(sprintf('Error in face areas at i=%d k=%d',i,k))
                %                end
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Calculate the signed distance from seed to face            
                % Find any point on face, for simplicity just choose the first vertex
                
                pt=Face_Vertices(1,:);
                
                % Get outward normal to face
                n=Face_Normals_i{j};
                
                % The signed distance to the face is
                D_spvals_val(r)=dot(pt-xi,n);
                
                % %%%%%%%%%%%%%%%%% For checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                point_on_face = sum(Face_Vertices)/length(Face_Vertices);
                
                % Calculate K(x,y) and K(x,z) for x = point on the face
                %                k_xk = kper_bx(point_on_face,xk,bx);
                
                % NOTE!
                % We think that k_xi should always be zero, possible speed up
                %                k_xi = kper_bx(point_on_face,xi,bx);
                
                % del_k = K(x,y) - K(x_z)
                % del_k keeps track of which periodic copies of the seeds are used to construct the face
            
                %                del_k = k_xk - k_xi;
                %                dist_ik = norm(xi-xk+del_k*diag(bx));
                
                %                Dk_spvals_val(r)=dist_ik;
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                % Increment r to keep track of number of non-zero entries of A and D
                r=r+1;
                
            end
        end
    end
    
    % Now use AD_spvals and A_spvals and D_spvals to construct the sparse matrices containing the face areas (A)
    % and the interparticle distances (D). A is symmetric, D is not as (i,j) measures the (signed) distance from
    % particle i to the face of cell i that is formed from the cell of a periodic copy of particle j. To get the
    % interparticle distance we need to form D+D'
    
    A=sparse(AD_spvals_i(1:r-1),AD_spvals_j(1:r-1),A_spvals_val(1:r-1),Nw,Nw);
    D=sparse(AD_spvals_i(1:r-1),AD_spvals_j(1:r-1),D_spvals_val(1:r-1),Nw,Nw);

    % Construct the symmetric interparticle distance matrix (works in periodic and non-periodic case)
    D=D+D';

    % Reciprocal distance (we have to do it this way because 1./D when D is sparse makes a full matrix with lots of NaNs
    inverse_D=spfun(@(x) 1./x, D);
    
    % Construct the off-diagonal entries
    H=-0.5*A.*inverse_D;
    
    % Now having all the off-diagonal entries we can add the diagonal entries
    H=H-spdiags(sum(H)',0,Nw,Nw);
  
    % We finished!
    success=true;
end


