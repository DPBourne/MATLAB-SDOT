function [w_0,exiterr,v_0,EXITFLAG] = SDOT_damped_Newton(w_0,X,target_vols,bx,periodic,percent_tol)
% [w_0,exiterr,v_0,EXITFLAG] = SDOT_damped_Newton(w_0,X,target_vols,bx,periodic,percent_tol)
%
% Note: Say what this functions does, i.e., solve the semi-discrete OT
% problem with cost ... , source measure ..., target measure ... Algorithm:
% damped Newton
%
% This function performs the damped Newton method of extremising g and
% returns the weights corresponding to a extremum
%
% Input arguments:
%
% w_0 is an Nx1 array containing the initial guess for the weights.
%
% X is an Nx3 array containing the x-, y- and z-coordinates of the seeds.
%
% target_vols is an Nx1 array containing the desired volumes of the cells.
%
% bx=[L1,L2,L3] is a 1-by-3 array containing the dimensions of a 
% rectangular box with vertices (0,0,0), (L1,0,0), (0,L2,0), (0,0,L3), 
% (L1,L2,0), (L1,0,L3), (0,L2,L3), (L1,L2,L3).
%
% periodic is a logical argument: If periodic=true, then the Laguerre
% diagram is periodic. If periodic=false, then the Laguerre diagram is 
% non-periodic.
%
% percent_tol is the tolerance for the convex optimization algorithm. This function
% produces a Laguerre diagram with cells of given volumes (target_vols) up 
% to tol percent error.
%
% Output arguments:
%
% w_0 is an Nx1 array of weights so that the (periodic or non-periodic)
% Laguerre diagram generated by (X,w) in the box of dimensions 
% (L1,L2,L3) centred at (L1/2,L2/2,L3/2) has Laguerre cells with volumes 
% target_vols up to tol percentage error.
%
% exiterr is the maximum percentage error of the volumes of the cells
%
% v_0 is a Nx1 array containing the actual volumes of the cells in 
% the Laguerre diagram generated by (X,w). These volumes agree with the
% desired volumes (target_vols) up to tol percentage error.
%
% EXITFLAG can take the value 0 or 1. 0 means that there is at least one 
% zero-volume cell. 1 means that the algorithm is successful.
%
% Last updated: 4 August 2021

[n,~] = size(X);

% To obtain a unique solution for the weights we impose the condition that w_0(n)=0
w_0(n) = 0;

% Volumes corresponding to initial weights, used to determine epsilon
[v_0,~,~] = mexPD(bx,X,w_0,periodic);

% Epsilon2 is defined as the minimum initial volume (corresponding to w_0)
epsilon2 = min(v_0);

% Test for zero volumes
if(epsilon2<=0)
    warning('With the w_0 specified, there is at least one zero-volume cell')
    EXITFLAG=0;
    w_0=[]; exiterr=[]; v_0=[];
    return
end

% The epsilon in Kitagawa, Merigot & Thibert (2019) is derived from the minimum of the target volumes
% and initial volumes (corresponding to w_0)
epsilon1 = min(target_vols);
epsilon = 0.5*min([epsilon1,epsilon2]);

% Given the percentage tolerance for the target volumes we calculate the absolute tolerance for the volumes
tol=epsilon1*percent_tol/100;

% Error using the current weights (v_0 in this case)
err_0=max(abs(v_0-target_vols));

% Preallocation (for speed)
v_k=zeros(n,1);

% Start the while loop
while(err_0>tol)
        
    % Obtain the gradient and Hessian
    [~,Dg,H,~]=kantorovich(w_0,X,target_vols,bx,periodic);
         
    % The Hessian H is singular, with one-dimensional kernel (1,1,...,1), truncate H to produce the non-singular H_mod
    H_mod = H(1:n-1,1:n-1);
    
    % Solve the linear system
    v_k(1:n-1) = -H_mod\Dg(1:n-1);
    v_k(n) = 0;
       
    % Backtracking
    backtracking=true;
    % l controls length of proposed Newton step, which is 2^(-l)*v_k
    l=0;
            
    while(backtracking)
        %disp(sprintf('Step %d',l));
        
        % Newton step: proposed new weights
        w_l=w_0+v_k*2^(-l);
        
        % Obtain the volumes corresponding to the proposed step
        [v_l,~,~]=mexPD(bx,X,w_l,periodic);
        
        % Find the minimum volume
        min_vol=min(v_l);
        
        % If the minimum volume is greater than epsilon we have a valid step (meaning not too close to a zero volume)
        if(min_vol>epsilon)
            
            % Check that there is a sufficient decrease in error
            % err_l is the error after the proposed step
            err_l=max(abs(v_l - target_vols));
            % Ratio of the error after proposed step compared to error before backtracking
            ratio=err_l/err_0;
            
            % Threshold
            threshold=1-2^(-(l+1));
            
            % If there is a sufficient decrease in error accept the Newton step with given l
            if(ratio<threshold)                
                % Accept the proposed weights as the new weights
                w_0=w_l;
                v_0=v_l;
                % Stop backtracking
                backtracking=false;

                % If the ratio > threshold then increment l (take a smaller Newton step)
            else
                l=l+1;
            end
            % If the minimum volume is too small (below epsilon) then take smaller Newton step increment l
        else
            l=l+1;
        end
    end
    % Update the current error
    err_0=max(abs(v_0-target_vols));
    
end

exiterr=max(abs(v_0-target_vols)./target_vols)*100;
EXITFLAG=1;

end