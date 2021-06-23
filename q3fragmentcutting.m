function rho = q3fragmentcutting(U1, V)
% fragment結合後の最終出力（密度演算子）

% Xの射影演算子
Xp1 = (1/sqrt(2))*[1;1];
Xn1 = (1/sqrt(2))*[1;-1];
Xp2 = Xp1*(Xp1)';
Xn2 = Xn1*(Xn1)';
X_list = [Xp2 Xn2];

%Yの射影演算子
Yp1 = (1/sqrt(2))*[1;1i];
Yn1 = (1/sqrt(2))*[1;-1i];
Yp2 = Yp1*(Yp1)';
Yn2 = Yn1*(Yn1)';
Y_list = [Yp2 Yn2];

%Zの射影演算子
Zp1 = [1;0];
Zn1 = [0;1];
Zp2 = Zp1*(Zp1)';
Zn2 = Zn1*(Zn1)';
Z_list = [Zp2 Zn2];

%Iの射影演算子
I_list = [Zp2 Zn2];

% rhoIを以下で別途計算
I = eye(2);
rho1I = zeros(2,2);
U1_0_U1 = U1 * [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] * (U1)';
for i=1:2
    Ii = kron(I,I_list(:,(2*i-1):2*i)); 
    Ii_U1_0_U1 = Ii * U1_0_U1;
    rho1I = rho1I + [Ii_U1_0_U1(1,1) + Ii_U1_0_U1(2,2)  Ii_U1_0_U1(1,3) + Ii_U1_0_U1(2,4); Ii_U1_0_U1(3,1) + Ii_U1_0_U1(4,2) Ii_U1_0_U1(3,3) + Ii_U1_0_U1(4,4)];
end

rho2I = zeros(4,4);
for j=1:2
    Ij = I_list(:,(2*j-1):2*j); 
    Ij_0 = kron(Ij,[1 0; 0 0]);
    V_Ij_0_V = V * Ij_0 * V';
    rho2I = rho2I + V_Ij_0_V;
end

rhoI = kron(rho1I, rho2I);

% ここまではrhoIの計算

rhoX = rhotoM(U1, V, X_list);
rhoY = rhotoM(U1, V, Y_list);
rhoZ = rhotoM(U1, V, Z_list);

rho = (rhoX + rhoY + rhoZ + rhoI)/2;
end