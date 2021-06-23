function rho1M = rho1toM(U1, M_list)
% rho1Mは最初のfragmentに関して、第2qubitをM（X, Y, Z）
% の基底(|0>, |+>など)に射影した後の第1qubitの密度演算子の出力
% U1は最初のfragmentの演算子（演算子の構成に注意）

I = eye(2);
rho1M = zeros(2,2);
U1_0_U1 = U1 * [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] * (U1)';
for i=1:2
    Mi = kron(I,M_list(:,(2*i-1):2*i)); 
    Mi_U1_0_U1 = Mi * U1_0_U1;
    rho1M = rho1M + (3-2*i) * [Mi_U1_0_U1(1,1) + Mi_U1_0_U1(2,2)  Mi_U1_0_U1(1,3) + Mi_U1_0_U1(2,4); Mi_U1_0_U1(3,1) + Mi_U1_0_U1(4,2) Mi_U1_0_U1(3,3) + Mi_U1_0_U1(4,4)];
end
end
