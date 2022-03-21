function yPrime = Mass_ty_Prime(t,y,yLag,varargin)
MCte      = varargin{:};
M         = MCte*(t+1)*(y(1)+2)*(y(2)-4);
yPrime(1) = y(2);
yPrime(2) = - yLag(1,1); 
yPrime    = M*yPrime(:);
yPrime    = yPrime(:);
