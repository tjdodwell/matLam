function [A,B,D,H] = makeABDH2(model)


% upper and lower coordinates
z = zeros(1,model.numPly+1);
z(1) = 0;
for i = 2:model.numPly+1
   z(i) = z(i-1) + model.t(i-1);
end
z = z - mean(z);

model.material.nu21=model.material.nu12*(model.material.E2/model.material.E1);
factor=1-model.material.nu12*model.material.nu21;

Q = zeros(5);
Q(1,1)=model.material.E1/factor;
Q(1,2)=model.material.nu12*model.material.E2/factor;
Q(2,1)=Q(1,2);
Q(2,2)=model.material.E2/factor;
Q(3,3)=model.material.G12;
Q(4,4)=model.material.SF*model.material.G23;
Q(5,5)=model.material.SF*model.material.G13;

%______________________________________________
A = zeros(5); B = zeros(5); D = zeros(5); H = zeros(5); T = zeros(5);
for k=1:model.numPly
    
    phi = model.ss(k);
    
    % Transformation Matrix    
    c = cos(phi); s = sin(phi);
    T(1,1) = c^2; T(1,2) = s^2; T(1,3) = 2*c*s;
    T(2,1) = s^2; T(2,2) = c^2; T(2,3) = -2*c*s;
    T(3,1) = -c*s; T(3,2) = c*s; T(3,3) = c^2-s^2;
    T(4,4) = c; T(4,5) = s;
    T(5,4) = -s; T(5,5) = c;
    
    % [Q] in structural axes
    invT = inv(T);
    Qbar= invT*Q*(invT');
    
    A= A + Qbar*(z(k+1)-z(k));
    B= B + Qbar*(z(k+1)^2-z(k)^2)/2;
    D= D + Qbar*(z(k+1)^3-z(k)^3)/3;
    H= H + Qbar*(z(k+1)-z(k));
    
end

A = A(1:3,1:3); B = B(1:3,1:3); D = D(1:3,1:3);
H = H(4:5,4:5);

end
