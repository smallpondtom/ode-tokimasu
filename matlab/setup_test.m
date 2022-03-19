addpath("exampleODE\");

switch which_test
    case 1
        ic = [0;0];
        tspan = [0,5];
        fs = @test01;
        fe = @exactSol01;
    case 2
        ic = [0.5];
        tspan = [0, 10];
        fs = @test02;
        fe = @exactSol02;
    case 3
        ic = [2;1];
        tspan = [0,20];
        fs = @test03;
        fe = @exactSol03;
    case 4
        ic = [0;2];
        tspan = [0,20];
        fs = @test04;
        fe = @exactSol04;
    case 5
        ic = [1;-8];
        tspan = [0,10];
        fs = @test05;
        fe = @exactSol05;
    case 6
        ic = [1];
        tspan = [1,12];
        fs = @test06;
        fe = @exactSol06;
    case 7
        ic = [3;0;0];
        tspan = [0,6*pi];
        fs = @test07;
        fe = @exactSol07;
    case 8
        ic = [0.4;0;0;2];
        tspan = [0,6*pi];
        fs = @test08;
        fe = @exactSol08;
    case 9
        ic = [-4];
        tspan = [1,10];
        fs = @test09;
        fe = @exactSol09;
    case 10
        ic = [1;1;1;1;2;0;1];
        tspan = [0,10];
        fs = @test10;
        fe = @exactSol10;
    case 11
        ic = [0.5;-1;1/3;-1/8;-35];
        tspan = [-1,12];
        fs = @test11;
        fe = @exactSol11;
    case 12
        epsilon = 1e-6;
        ic = [-0.1/(epsilon + 0.01)^0.5,epsilon / (epsilon + 0.01)^0.5];
        tspan = [-0.1,0.1];
        fs = @test12;
        fe = @exactSol12;
    otherwise
        error("Select only tests 1 to 12.");
end