function [g,Dg,H,actual_vols] = get_g_kper(w,X,target_vols,bx)

% [g,Dg,H,actual_vols] = get_g(w,X,target_vols,bx)
%
% WARNING! Periodic diagrams only
%
% This function computes g, Dg, H for periodic Laguerre diagrams using the
% function kper_bx
%
% Input arguments
%         w            - is the weights, an Nx1 array
%         X            - is the positions, an Nx3 array
%         target_vols  - target volumes, an Nx1 array
%         bx           - box size, a 1x3 array
%
% Return arguments
%         g            - function g(w;x)
%         Dg           - gradient of g wrt w
%         H            - the Hessian d^2g, a sparse matrix
%         actual_vols  - the actual volumes of the Laguerre diagram with seeds X and weights w
    
    %% Catch errors
    
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
    % Note the diagram is assumed periodic so the last argument must be true
    [actual_vols,transport_costs,~,vfn]=mexPDallfaces(bx,X,w,true);

    % Gradient of g
    Dg=actual_vols-target_vols;

    % Definition of g, a convex function of the weights 
    g=dot(Dg,w)-sum(transport_costs); 
        
    % FUTURE DEVELOPMENT? Preallocate memory for the sparse matrix entries
    % INSERT HERE
    
    % Matrices to store the values and indices for the sparse Hessian matrix
    H_spvals_i=[];
    H_spvals_j=[];
    H_spvals_val=[];
    
    % Index to keep track of the number of entries in the sparse matrix
    r=1;
    
    % we now start computing the hessian values
    for i=1:Nw
        % Coordinates of the seed location of the ith cell
        xi=X(i,:);
        
        % The list of neighbour indices of the ith cell - boundaries have negative indices
        N_i=vfn{i,3};
        
        % The vertex indices for each of the faces between the cell and the neighbours in N_i
        Face_Index_i=vfn{i,2}(:,1);
        Face_Areas_i=vfn{i,2}(:,2);
        
        % The vertices for the ith cell
        Verticies_Cell_i=vfn{i,1};
        
        % The number of neighbours
        Num_N_i = size(N_i);
        
        % calculate H_ij for j in N_i 
        for j=1:Num_N_i
            
            % The cell index of neighbour j
            k = N_i(j);
            
            % We only need to calculate if it is a true neighbour (and we also exploit that the matrix is symmetric)
            
            if (k>i)  % Upper triangular elements only
                
                % The indices of the vertices for the face between cell i and cell k (neighbour j)
                Face_Vertex_Indices = Face_Index_i{j};
                
                % The face vertices for the face between cell i and cell k (neighbour j)
                Face_Vertices = Verticies_Cell_i(Face_Vertex_Indices,:);
                                           
                % Coordinates of the the seed for cell k
                xk = X(k,:);
                                
                % Find any point in the interior of the face
                point_on_face = sum(Face_Vertices)/length(Face_Vertices);
                
                % Calculate K(x,y) and K(x,z) for x = point on the face
                k_xk = kper_bx(point_on_face,xk,bx);
                
                % NOTE!
                % We think that k_xi should always be zero, possible speed up
                k_xi = kper_bx(point_on_face,xi,bx);
                
                % del_k = K(x,y) - K(x_z)
                % del_k keeps track of which periodic copies of the seeds are used to construct the face
            
                del_k = k_xk - k_xi;
                
                % Area of the face
                % OLD METHOD, uses MatGeom library, very slow
                % a_ik_polyarea = polygonArea3d(Face_Vertices);
              
                a_ik = Face_Areas_i{j};
                
                % The correct distance between neighbours accounts for the fact that the neighbour could
                % be a periodic copy of the original seed
                
                dist_ik = norm(xi-xk+del_k*diag(bx));
                
                % Indices and value of the Hessian
                H_spvals_i(r)=i;
                H_spvals_j(r)=k;
                H_spvals_val(r)=-0.5*a_ik/dist_ik;
                r=r+1;
            end            
        end
        
    end
    
    % Now use H_spvals to construct the sparse matrix
    H=sparse(H_spvals_i,H_spvals_j,H_spvals_val,Nw,Nw);
    
    % Having calculated the upper triangular part we can find all off-diagonal entries by adding the transpose
    H=H+H';

    % Now having all the off-diagonal entries we can add the diagonal entries
    H=H-spdiags(sum(H)',0,Nw,Nw);
    
end
