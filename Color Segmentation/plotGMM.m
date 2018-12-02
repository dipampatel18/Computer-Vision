% Reference: https://www.mathworks.com/matlabcentral/fileexchange/16543-plot_gaussian_ellipsoid
function plotGMM(mu, covar)
[x,y,z] = sphere(20);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(covar); 
if any(d(:) < 0)
   d = max(d,0);
end
d = 1 * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(mu, 1, size(ap,2)); 
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));
surf(gca, xp,yp,zp);
end