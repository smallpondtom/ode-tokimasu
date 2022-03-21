function yPrime = MassCte_Prime(t,y,yLag,varargin)
MCte      = varargin{:};
yPrime(1) = y(2);
yPrime(2) = - yLag(1,1); 
yPrime    = (MCte)*yPrime(:);
yPrime    = yPrime(:);
