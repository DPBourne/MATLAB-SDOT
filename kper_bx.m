function [k,periodic_distance] = kper_bx(x,y,bx)
%
% function [k,periodic_distance] = kper_bx(x,y,bx)
%
% This function computes the periodic distance between x and y, where the
% periodic distance between two points x and y in R^d is defined as
%
% periodic_distance = \min_{l \in \mathbb{Z}^d} |x-y+Bl|
%
% and
%
% k = argmin_{l \in \mathbb{Z}^d} |x-y+Bl|
%    
% where B=diag(bx).
%
% This function returns k and periodic_distance
%
% Inputs are x  - a 1xN vector
%            y  - a 1xN vector
%            bx - a 1xN vector of box side lengths
%
% Outputs    k                 - a 1xN vector of integers
%            periodic_distance - a double, the periodic distance between x and y in a periodic box size given by bx

% Calculate the dimension d
[~,d]=size(x);

% Diagonal matrix, bx (box sizes) on the diagonal
B=diag(bx);Binv=diag(1./bx);

% Calculate k
k=floor((y-x)*Binv+0.5*ones(1,d));
periodic_distance=norm(x-y+k*B);
    
end 