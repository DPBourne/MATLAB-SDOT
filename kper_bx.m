function [k] = kper_bx(x,y,bx)
% This caculates the K in the defintion of periodic distance,

% Inputs are x  - a 1xN vector
%            y  - a 1xN vector
%            bx - a 1xN vector of box side lengths

% Outputs k a 1xN vector of integers
% $$$     
     [mx,nx]=size(x);
% $$$     if(mx~=1)
% $$$         error('x should be a 1xN vector');
% $$$     end
% $$$     
     [my,ny]=size(y);
% $$$     if(my~=1)
% $$$         error('y should be a 1xN vector');
% $$$     end
% $$$    
% $$$     if(ny~=nx)
% $$$         error('x and y should both be 1xN vectors');
% $$$     end
% $$$     
     [mbx,nbx]=size(bx);
% $$$     if(mbx~=1)
% $$$         error('bx should be a 1XN vector');
% $$$     end
% $$$     
% $$$     if(nbx~=nx)
% $$$         error('bx,x and y should be 1xN vectors');
% $$$     end
    
    % Diagonal matrix made of bx (box sizes)
    Binv=diag(1./bx);
    k = floor((y-x)*Binv+0.5*ones(1,nx));
    
end 