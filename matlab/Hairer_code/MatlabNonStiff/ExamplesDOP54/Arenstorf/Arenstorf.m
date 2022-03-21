function dx = Arenstorf(t,x,Param)
% See Hairer université de Genève
% Driver for DOPRI5
% SUBROUTINE FAREN(N,X,Y,F,RPAR,IPAR)
% C --- ARENSTORF ORBIT
%         IMPLICIT REAL*8 (A-H,O-Z)
%         DIMENSION Y(N),F(N),RPAR(2)
%         AMU=RPAR(1)
%         AMUP=RPAR(2)
%         F(1)=Y(3)
%         F(2)=Y(4)
%         R1=(Y(1)+AMU)**2+Y(2)**2
%         R1=R1*SQRT(R1)
%         R2=(Y(1)-AMUP)**2+Y(2)**2
%         R2=R2*SQRT(R2)
%         F(3)=Y(1)+2*Y(4)-AMUP*(Y(1)+AMU)/R1-AMU*(Y(1)-AMUP)/R2
%         F(4)=Y(2)-2*Y(3)-AMUP*Y(2)/R1-AMU*Y(2)/R2
%         RETURN
%         END 
% RPAR(1) = 0.012277471D0
% RPAR(2) = 1.D0-RPAR(1)
% Y(1)=0.994D0
% Y(2)=0.0D0
% Y(3)=0.0D0
% Y(4)=-2.00158510637908252240537862224D0
% X=0.0D0
% XEND=17.0652165601579625588917206249D0
% RTOL=1.0D-7
% ATOL=RTOL

% ------------
Amu  = Param(1);
Amup = Param(2);
dx    = zeros(4,1);
dx(1) = x(3);
dx(2) = x(4);
R1    = (x(1)+Amu)^2 + x(2)^2;
R1    = R1*sqrt(R1);
R2    = (x(1)-Amup)^2 + x(2)^2;
R2    = R2*sqrt(R2);
dx(3) = x(1) + 2*x(4) - Amup*(x(1)+Amu)/R1 - Amu*(x(1)-Amup)/R2;
dx(4) = x(2) - 2*x(3) - Amup*x(2)/R1 - Amu*x(2)/R2;



