function K = computeGlobalStiffness(model,msh)

%  computeGlobalStiffness.m ~ Dr. Tim Dodwell

%  Construct global stiffness matrices K for solution of linear system of equations
%                    K * d = f
%  where d = (w_1, \phi_x1, \phi_y1, ..., w_N, \phi_xN, \phi_yN)'.
%
%% Construct ABD and H matrices for each element
%  Element-wise construction of plate thickness, ply angle interpolation,
%  and constitutive relation matrices

[A, B,  D,  H] = makeABDH2(model);

%% Allocate storage for global stiffness matrices

indx_j  = repmat(1:msh.nedof,msh.nedof,1); 
indx_i  = indx_j';
Kindx.i = msh.e2g(:,indx_i(:)); 
Kindx.j = msh.e2g(:,indx_j(:));

Ke_all  = zeros(msh.nedof^2,msh.nel); 

%% Construct stiffness matrices for transverse calculation

%  Loop over elements:

for ie = 1:msh.nel
    
    K_elem = zeros(msh.nedof);
    
    % Loop over integration points
    
    for ip = 1:4
        
        % Load shape functions and their derivatives (wrt global coordinate
        % system) for integration point ip, as well as Jacobian determinant
        
        [~,dNdX,detJ] = elementShapeFunctions(msh,1,ip,'full');
        
        % Construct [B_b], matrix bending (called G_1 in notes):
        
        B_b = zeros(msh.dof,msh.nedof);
        B_b(1,msh.nnodel+1:2*msh.nnodel)   = dNdX(:,1);
        B_b(2,2*msh.nnodel+1:3*msh.nnodel) = dNdX(:,2);
        B_b(3,msh.nnodel+1:2*msh.nnodel)   = dNdX(:,2);
        B_b(3,2*msh.nnodel+1:3*msh.nnodel) = dNdX(:,1);
        
        % Compute element stiffness matrices:
        % Deff ~ decoupling membrane and bending problems. Take D* = D - B'inv(A)B (Schur complement)
        
        Deff = D - (B' * inv(A) * B);
        
        % Construct stiffness matrix:
        % (Gauss integration weight is equal to one)
        
        K_elem = K_elem + (B_b' * Deff * B_b) * detJ; 
        
    end
    
    G_elem = zeros(msh.nedof);
   
    ip     = 1;
    
    % Load shape functions and their derivatives (wrt global coordinate
    % system) for integration point ip, as well as Jacobian determinant
    
    [Ni,dNdX,detJ] = elementShapeFunctions(msh,ie,ip,'reduced');
  
    % Construct [B_s], matrix shear (called G_2 in notes):
    
    B_s = zeros(2,msh.nedof);
    B_s(1,1:msh.nnodel)                 = dNdX(:,1);
    B_s(2,1:msh.nnodel)                 = dNdX(:,2);
    B_s(1,msh.nnodel+1:2*msh.nnodel)    = Ni;
    B_s(2,2*msh.nnodel+1:3*msh.nnodel)  = Ni;
    
    % Compute element stiffness matrices:
   
    K_elem = K_elem + (B_s'* H *B_s) * msh.IP_w(ip) * detJ;
   
    % Store element stiffness matrices as columns
    
    Ke_all(:, ie) = K_elem(:);
    
end

%% Write data into global storage
K = sparse(Kindx.i',Kindx.j',Ke_all);

end



