function Mass = MassFcn_ty(t,y)
Mass = [2,-3;-1,4]*(t+1)*(y(1)+2)*(y(2)-4);