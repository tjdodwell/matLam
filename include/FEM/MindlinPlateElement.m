function Ke = MindlinPlateElement(ie,msh,mat)

% Mindlin Plate - with selective integration to reduce shear locking

% Degrees of freedom per element d = [u,v,w,tx,ty] - 5 dofs per node * 4 per element = 20

K_elem = zeros(20); % Initialise Element Stiffness Matrix

ne = msh.nnodel; % store nodes per element (4) to save writing out all the time

% Loop over integration points

for ip = 1:4
        
        % Load shape functions and their derivatives (wrt global coordinate
        % system) for integration point ip, as well as Jacobian determinant
        
        [~,dNdX,detJ] = elementShapeFunctions(msh,1,ip,'full');

        % Construct [B], in plane matrix

        B = zeros(6,msh.nedof);

        B(1,1 : ne) = dNdX(:,1); % du/dx
        B(2,ne + 1 : 2 * ne) = dNdX(:,2); % dv/dx
        B(3,1 : ne) = dNdX(:,2);
        B(3,ne + 1 : 2 * ne) = dNdX(:,1);

        B(4,3*ne + 1 : 4*ne) = dNdX(:,1); % dtx / dx
        B(5,4*ne + 1 : 5*ne) = dNdX(:,2); % dty / dy
        
        B(6,3*ne + 1 : 4*ne) = dNdX(:,2); % dtx / dy
        B(6,4*ne + 1 : 5*ne) = dNdX(:,1); % dty / dx
        
        % Compute element stiffness matrices:
      
        % Construct stiffness matrix:
        % (Gauss integration weight is equal to one)

        C = [mat.A, mat.B; mat.B', mat.D];

        
        K_elem = K_elem + (B' * C * B) * detJ; % Assumes weight is 1
        
end
   
    % Reduced integration on shear terms to prevent shear locking

    ip  = 1; 
    
    % Load shape functions and their derivatives (wrt global coordinate
    % system) for integration point ip, as well as Jacobian determinant
    
    [Ni,dNdX,detJ] = elementShapeFunctions(msh,ie,ip,'reduced');
  
    % Construct [B_s], matrix shear:
    
    B_s = zeros(2,msh.nedof);
    B_s(1,2*ne + 1 : 3*ne)                 = dNdX(:,1); % dw/dx
    B_s(2,2*ne + 1 : 3*ne)                 = dNdX(:,2); % dw/dy
    B_s(1,3*ne + 1 : 4*ne)                 = Ni;
    B_s(2,4*ne + 1 : 5*ne)                 = Ni;
    
    % Compute element stiffness matrices:
   
    K_elem = K_elem + (B_s'* mat.H *B_s) * 4.0 * detJ; % Reduced integration so 4 is integration weight
   
    % Store element stiffness matrices as columns

    Ke = K_elem(:);

end