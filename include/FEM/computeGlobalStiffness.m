function K = computeGlobalStiffness(model,msh)

%  computeGlobalStiffness.m ~ Dr. Tim Dodwell

%  Construct global stiffness matrices K for solution of linear system of equations
%                    K * d = f
%  where d = (w_1, \phi_x1, \phi_y1, ..., w_N, \phi_xN, \phi_yN)'.
%
%% Construct ABD and H matrices for each element
%  Element-wise construction of plate thickness, ply angle interpolation,
%  and constitutive relation matrices

mat = makeABDH2(model);

%% Allocate storage for global stiffness matrices

indx_j  = repmat(1:msh.nedof,msh.nedof,1); 
indx_i  = indx_j';
Kindx.i = msh.e2g(:,indx_i(:)); 
Kindx.j = msh.e2g(:,indx_j(:));

Ke_all  = zeros(msh.nedof^2,msh.nel); 

%% Construct stiffness matrices for transverse calculation

%  Loop over elements:

switch lower(model.type)

    case 'mindlin'
        for ie = 1:msh.nel
            Ke_all(:, ie) = MindlinPlateElement(ie,msh,mat);
        end

    case 'zigzag'
        for ie = 1:msh.nel
            Ke_all(:, ie) = ZigZagPlateElement(ie,msh,mat);
        end

end

%% Write data into global storage
K = sparse(Kindx.i',Kindx.j',Ke_all);

end



