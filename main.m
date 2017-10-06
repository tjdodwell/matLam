close all
clear all

addpath('include/FEM/');
addpath('include/preProcessing/');
addpath('include/postProcessing/');

model = modelSetup();
msh   = makeMesh(model);

K = computeGlobalStiffness(model,msh);

% Lets do boundary conditions

msh.lhs = find(msh.coords(:,1) < 1e-6);
msh.rhs = find(msh.coords(:,1) > model.Lx - 1e-6);

bnd_left  = [msh.lhs; msh.lhs + msh.nnod; msh.lhs + 2 * msh.nnod];
bnd_right = [msh.rhs; msh.rhs + msh.nnod; msh.rhs + 2 * msh.nnod];

bnd = [bnd_left; bnd_right];

free = 1 : msh.tdof; free(bnd) = [];

% Solve the problem

U1 = zeros(msh.tdof,1); U0 = zeros(msh.tdof,1);

t = 0.0;

for i = 1 : model.timesteps
    
    t = t + model.dt;

    U1(msh.rhs) = model.A * sin(model.omega * t);

    U1(free) = K(free,free) \ (-K(free,bnd) * U1(bnd));

    % Plot Solution

    scalar_point.name = 'displacement';
    scalar_point.data = U1(1:msh.nnod);

    matlab2vtk(strcat('solution_',int2str(i),'.vtk'),'DMTA', msh, 'quad', [], scalar_point, []);

end

