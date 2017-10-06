function [Ni,dNdX,detJ] = elementShapeFunctions(msh,ie,ip,integration_option)
    
    switch lower(integration_option);
        
        case 'full'

         [IP_X,IP_W] = ip_quad;
    
     [N, dNdu] = shapeFunctionQ4(IP_X);

        Ni = N{ip}; dNdui = dNdu{ip};
    
        case 'reduced'
            
            Ni = msh.N{1}; dNdui = msh.dNdu{1};
         
     end

    J           = msh.coords(msh.elements(ie,:),:)'*dNdui';
    detJ        = det(J);
                
    dNdX        = dNdui'*inv(J);
    
end

function [N, dNdu] = shapeFunctionQ4(IP_X)
% TJD - June 2014
nip  = 4;
N    = cell(nip,1);
dNdu = cell(nip,1);
    for i = 1:nip
        xi = IP_X(i,1); eta = IP_X(i,2);
        shp=0.25*[ (1-xi)*(1-eta);
                   (1+xi)*(1-eta);
                   (1+xi)*(1+eta);
                   (1-xi)*(1+eta)];
        deriv=0.25*[-(1-eta), -(1-xi);
                     1-eta,    -(1+xi);
                     1+eta,      1+xi;
                     -(1+eta),   1-xi];
        N{i} = shp;
        dNdu{i} = deriv';
    end
end % end function shapeFunctionQ4

function [IP_X,IP_W] = ip_quad
% TJD - June 2014
% Gauss quadrature for Q4 elements
% option 'complete' (2x2)
% option 'reduced'  (1x1)
% nip: Number of Integration Points
% ipx: Gauss point locations
% ipw: Gauss point weights
            IP_X=...
              [ -0.577350269189626 -0.577350269189626;
                 0.577350269189626 -0.577350269189626;
                 0.577350269189626  0.577350269189626;
                -0.577350269189626  0.577350269189626];
            IP_W=[ 1;1;1;1]; 
end % end of function ip_quad

