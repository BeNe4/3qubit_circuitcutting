function rhoM = rhotoM (U1, V, M_list)
% rhoMはfragment結合後のMに関する密度演算子
% e.g. rhoX はrho1toXとrho2toMのテンソル積
% 演算子の順番に注意（U1の構成に関して、演算の順番はCX,HIなのかHI,CXなのか、等...）

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

rho1M = rho1toM(U1, M_list);
rho2M = rho2toM(V, M_list);
rhoM = kron(rho1M, rho2M);

end