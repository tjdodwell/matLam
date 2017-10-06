function [A, B, D, H] = constructZigZagMatrices(model)


% Compute the interfaces

% upper and lower coordinates
z = zeros(1,model.numPly+1);
z(1) = 0;
for i = 2:model.numPly+1
   z(i) = z(i-1) + model.t(i-1);
end
z = z - mean(z);

% Compute Q - composite matrix in local axis

model.material.nu21=model.material.nu12*(model.material.E2/model.material.E1);
factor=1-model.material.nu12*model.material.nu21;

Q = zeros(5);
Q(1,1)=model.material.E1/factor;
Q(1,2)=model.material.nu12*model.material.E2/factor;
Q(2,1)=Q(1,2);
Q(2,2)=model.material.E2/factor;
Q(3,3)=model.material.G12;

G23 = model.material.G23;
G13 = model.material.G13;

G13i = model.material.G13i;
G23i = model.material.G23i;

E_int = model.material.E2;
nu12_int = model.material.nu12;


% Compute Zig-Zag Matrices

[G1,	G2] = computeGs(model,G13,G23,G13i,G23i); % Compute laminate shear moduli

A = zeros(3); B = zeros(3,7); D = zeros(7); H = zeros(4);

for k = 1 : model.numPly

	phi = model.ss(k); 

	if (phi > 0) % This is a composite ply

		% Transformation Matrix    
    	c = cos(phi); s = sin(phi);
    	T = zeros(3);
	    T(1,1) = c^2; T(1,2) = s^2; T(1,3) = 2*c*s;
	    T(2,1) = s^2; T(2,2) = c^2; T(2,3) = -2*c*s;
	    T(3,1) = -c*s; T(3,2) = c*s; T(3,3) = c^2-s^2;
    
	    % Rotate [Q] to structural axes
	    invT = inv(T);
	    Qk = invT * Q * (invT');

	    % Shear Matrix

	    Hk = [cos(phi)^2 * G23 + sin(phi)^2 * G13, sin(phi) * cos(phi) * (G13 - G23);
			  sin(phi) * cos(phi) * (G13 - G23), (cos(phi) ^2) * G13 + (sin(phi) ^ 2) * G23];

	else % This is an interface layer

		% In this case notation is required.

		factor = 1.0 - nu12_int * nu12_int;

		Qk = [E_int/factor, nu12_int * E_int/factor, 0.0; nu12_int * E_int/factor, E_int/factor, 0.0; 0.0, 0.0, E_int / (2.0 * (1.0 + nu12_int))];

		Hk = [G23i,0.0;0.0,G13i];

	end

	% Compute A matrix - constant within each layer

	A = A + Qk * model.t(k);

	% Compute B matrix - linear in each layer - compute exactly with trapezoidal rule

	Bk0 = calB_phi(z(k),k,G13,G12,G13i,G23i,model);
	Bk1 = calB_phi(z(k+1),k,G13,G12,G13i,G23i,model)

	B = B + 0.5 * model.t(k) * Qk * (Bk1 + Bk0);

	% Compute D matrix - since quadratic in each layer - compute exactly with simpsons rule

	Bkhalf = calB_phi(0.5*(z(k)+z(k+1)),k,G13,G12,G13i,G23i,model);

	D = D + (1.0 / 6.0) * model.t(k) * (Bk0' * Qk * Bk0 + 4.0 * Bkhalf' * Qk * Bkhalf + Bk1' * Qk * Bk1);

	% Compute G matrix

	[beta1,	beta2] = computeBetas(k,G1,G2,G13,G23,G13i,G23i,model);

	Bb = [1.0, beta2, 0.0, 0.0; 0.0, 0.0, 1.0, beta1];

	G = G + model.t(k) * (Bb' * Qk * Bb) * model.t(k);

end % end for each ply

end



% model.t - contains layer thickness
% model.ss - contains stacking sequence


% Need a function which calculates Bphi

function [G1,	G2] = computeGs(model,G13,G23,G13i,G23i)
	G1 = 0.0; G1 = 0.0;
	for k = 1 : 2 : model.numPly
		phi = model.ss(k);
		Q11k = (cos(phi) ^ 2) * G13 + (sin(phi) ^ 2) * G23;
		Q22k = (cos(phi) ^ 2) * G23 + (sin(phi) ^ 2) * G13;
		G1 = G1 + model.t(k) / Q11k;
		G2 = G2 + model.t(k) / Q22k
	end
	for k = 2 : 2 : model.numPly % For the interfaces
		G1 = G1 + model.t(k) / G11i
		G2 = G2 + model.t(k) / G22i
	end
	G1 = G1 / sum(model.t); G2 = G2 / sum(model.t);
	G1 = 1 / G1; G2 = 1 / G2;
end

function [beta1,	beta2] = computeBetas(k,G1,G2,G13,G23,G13i,G23i,model)
	if (model.ss(k) < 0.0) % It is an interface
		beta1 = G1 / G13i - 1.0;
		beta2 = G2 / G23i - 1.0;
	else
		phi = model.ss(k);
		Q11k = (cos(phi) ^ 2) * G13 + (sin(phi) ^ 2) * G23;
		Q22k = (cos(phi) ^ 2) * G23 + (sin(phi) ^ 2) * G13;
		beta1 = G1 / Q11k - 1.0;
		beta2 = G2 / Q22k - 1.0;
	end
end


function B = calB_phi(z,k,G13,G12,G13i,G23i,model)
	B = zeros(3,7);
	[phi1,	phi2] = calPhi(z,k,G13,G12,G13i,G23i,model)
	B(1,1) = z; B(1,2) = phi1;
	B(2,3) = z; B(2,4) = phi2;
	B(3,5) = z; B(3,6) = phi1; B(3,7) = phi2;
end

function [phi1,	phi2] = calPhi(z,k,G13,G12,G13i,G23i,model)
	% Note that this is a linear function in z
	h = 0.5 * sum(model.t);
	phi = model.ss(k);
	Q11k = (cos(phi) ^ 2) * G13 + (sin(phi) ^ 2) * G23;
	Q22k = (cos(phi) ^ 2) * G23 + (sin(phi) ^ 2) * G13;
	phi1 = (z + h) * (G1 / Q11k - 1.0)
	phi2 = (z + h) * (G2 / Q22k - 1.0)
	if (k > 0)
		for i = 2 : model.numPly
			phi = model.ss(i);
			if(phi < 0.0)
				Q11i = G13i;
				Q22i = G23i;
			else		
				Q11i = (cos(phi) ^ 2 * G13) + (sin(phi) ^ 2) * G23;
				Q22i = (cos(phi) ^ 2 * G23) + (sin(phi) ^ 2) * G13;
			end
			phi1 = phi1 +  model.t(i-1) * (G1 / Q11i - G1 / Q11k)
			phi2 = phi2 +  model.t(i-1) * (G2 / Q22i - G2 / Q22k)
		end
	end
end



