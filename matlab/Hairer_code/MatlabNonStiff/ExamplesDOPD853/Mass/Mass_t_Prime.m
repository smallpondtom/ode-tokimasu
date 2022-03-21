function yPrime = Mass_t_Prime(t,y,yLag,varargin)
MCte      = varargin{:};
M         = MCte*(t+1);
yPrime(1) = y(2);
yPrime(2) = - yLag(1,1); 
yPrime    = M*yPrime(:);
yPrime    = yPrime(:);
