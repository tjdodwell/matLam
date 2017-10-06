function msh = makeMesh(model)
% -----------------------------------------------------------------------
% This code is released under GNU LESSER GENERAL PUBLIC LICENSE v3 (LGPL)
%
% Details are provided in license.txt file in the main directory
%
% 1/8/14 - Dr T. J. Dodwell - University of Bath - tjd20@bath.ac.uk
% -----------------------------------------------------------------------

% makemsh.m - Written (TJD - 3/6/2014)
%
% Creates Coarse Quadrilateral msh on [0,Lx] by [0,Ly] and refines uniformly to desired msh size
%

% --------------------------------
% (1)	Set up Coarse Rectangle
% --------------------------------
msh.coords = [0	0;
			   model.Lx 0;
			   model.Lx model.Ly;
			   0 model.Ly];
msh.elements = 1:4;
%
nodesOfRefinement = zeros(9,1);

edgeTable = zeros(0,3);  % Create an empty matrix with two columns


for ii = 1:model.meshRefinement % For each refinement
	
	visitedEdges = 0;
    nelem = 0;
	inode = 0;
    newcoords = [];
    nodesPreviousRefinement = size(msh.coords(:,1),1);
   
	for ie = 1:size(msh.elements,1); % Each Element

		nodesOfRefinement(1:4) = msh.elements(ie,:);

		% First Add Mid Point
		inode=inode+1;
		newcoords(inode,:) = 0.5*(msh.coords(msh.elements(ie,1),:) + msh.coords(msh.elements(ie,3),:));
		nodesOfRefinement(5)=inode + nodesPreviousRefinement;

		for edge = 1:4 % For Each Edge

			n1 = msh.elements(ie,edge); n2 = msh.elements(ie,mod(edge,4)+1);

			% Has Edge been visited before - if so return id = 1 and the node
			[id,oldNode] = edgeVisited(edgeTable,n1,n2);
			
			if id == 0 % If new edge add midpoint as new node
				inode=inode+1;
                nodesOfRefinement(5+edge) = inode + nodesPreviousRefinement;
				newcoords(inode,:) = 0.5*(msh.coords(n1,:)+msh.coords(n2,:));
				edgeTable(visitedEdges+1,1:2) = [n1,n2];
				edgeTable(visitedEdges+1,3) = inode + nodesPreviousRefinement;
				visitedEdges=visitedEdges+1;
			else
				nodesOfRefinement(5+edge) = oldNode;
			end

		end % For each edge

		newelements(nelem+1,:)=nodesOfRefinement([1,6,5,9]);
		newelements(nelem+2,:)=nodesOfRefinement([6,2,7,5]);
		newelements(nelem+3,:)=nodesOfRefinement([5,7,3,8]);
		newelements(nelem+4,:)=nodesOfRefinement([9,5,8,4]);
		nelem=nelem+4;

	end

	msh.elements = newelements;
	msh.coords = [msh.coords;newcoords];

end % for each refinelement

% Mesh Constructed

msh.nnod        = size(msh.coords,1);
msh.nel         = size(msh.elements,1);

msh.ndim        = 2;

if (strcmp('zigzag_reduced',model.type))
	msh.dof = 5;
else
	msh.dof = 3; % Default is Timoshenko
end



msh.dof         = 3;
msh.nnodel      = 4; % Nodes per Element
msh.nedof       = msh.nnodel*msh.dof;
msh.tdof        = msh.dof*msh.nnod;


msh.e2g = zeros(msh.nel,msh.nedof);

for ie = 1:msh.nel % For each element
		msh.e2g(ie,:) = [msh.elements(ie,:),msh.nnod + msh.elements(ie,:),2*msh.nnod + msh.elements(ie,:)];
end

[msh.nip,msh.IP_X, msh.IP_w]    = ip_quad(model.integrationOption);                   
[msh.N, msh.dNdu]        = shapeFunctionQ4(msh.IP_X);


end

function [id,node] = edgeVisited(edgeTable,n1,n2)

id = 0; node = 0;

numEdgesVisited = size(edgeTable,1);

temp = find(edgeTable(:,1) == n1);
for i = 1:length(temp)
	if edgeTable(temp(i),2) == n2
		id = 1;
		node = edgeTable(temp(i),3);
	end
end
if id == 0
	temp = find(edgeTable(:,2) == n1);
	for i = 1:length(temp)
		if edgeTable(temp(i),1) == n2
			id = 1;
			node = edgeTable(temp(i),3);
		end
	end
end

end

function [N, dNdu] = shapeFunctionQ4(IP_X,nnodel)
% TJD - June 2014
nip  = size(IP_X,1);
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

function [nip,IP_X,IP_W] = ip_quad(option)
% TJD - June 2014
% Gauss quadrature for Q4 elements
% option 'complete' (2x2)
% option 'reduced'  (1x1)
% nip: Number of Integration Points
% ipx: Gauss point locations
% ipw: Gauss point weights
    switch option
        case 'complete'
        nip = 4;
            IP_X=...
              [ -0.577350269189626 -0.577350269189626;
                 0.577350269189626 -0.577350269189626;
                 0.577350269189626  0.577350269189626;
                -0.577350269189626  0.577350269189626];
            IP_W=[ 1;1;1;1]; 
        case 'reduced'
            nip = 1;
            IP_X=[0 0];
            IP_W=[4];
    end % end of switch 'option'
end % end of function ip_quad

