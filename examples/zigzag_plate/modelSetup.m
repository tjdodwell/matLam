function model = modelSetup

model.plotme = true;
model.vtk_filename_base = 'results/solution_';

model.type = 'zigzag';

model.Lx = 10.0; % Length (x dimension)
model.Ly = 5.0; % Width (y dimension)



model.numPly = 1; % number of plys
model.t = [0.2];
model.ss	 = [0.0];

model.material.G13i = 2.0;
model.material.G23i = 2.0;

model.material.E1 = 130; % kN/mm^2
model.material.E2 = 9.25; % kN/mm^2
model.material.G12 = 5.1; % kN/mm^2
model.material.G23 = 5.13; % kN/mm^2
model.material.G13 = 5.13; % kN/mm^2
model.material.nu12 = 0.36;
model.material.SF = 5/6; % Shear Correction Factor
model.material.density = 1584E-9; % probably consistent units re: the above

model.timesteps = 10;
model.dt = 0.1;
model.A = 1.0;
model.omega = 1.0;

model.meshRefinement = 4;

model.integrationOption = 'reduced'; % Type of Integration, if you use 'complete' you get shear locking

end