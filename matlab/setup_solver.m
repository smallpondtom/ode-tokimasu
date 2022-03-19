addpath("unmodified\");
addpath("modified\");

switch which_method
    case 1
        solver = @Euler;
    case 2
        solver = @RK4;
    case 3
        solver = @RKF45;
    case 4
        solver = @DOPRI54;
    otherwise
        error("Select only tests 1 to 4.");
end