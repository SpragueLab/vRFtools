% vRF_2dGaussian_nonlinear
%
% function used to generate 2D nonlinear Gausisan 'filter' quickly &
% efficiently 
%
% INPUT:
% - eval_at: points (x,y) for evaluating the function at
% - params:  parameters for the function (mu_x, mu_y, sigma, exp), as
%            columns
%
% OUTPUT:
% - img: each row is a filter evaluated at eval_pts, for the corresponding
%        row of params' parameters
%
% NOTE: exp param is irrelevant here...that only comes into play after
% combining w/ stimulus mask...

function img = vRF_2dGaussian_nonlinear(eval_at,params)

% for efficiency, we'll run w/ eval_at as a row vector
if size(eval_at,2)==2
    eval_at = eval_at.'; % transpose if submitted as column vec
end

tmpr = hypot(eval_at(1,:) - params(:,1),eval_at(2,:) - params(:,2));

% exp(  -r^2 / (2*sigma^2)  )
img = exp( -tmpr.^2 ./ (2*params(:,3).^2) );



return;