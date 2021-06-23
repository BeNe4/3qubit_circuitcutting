function rho2M = rho2toM (V, M_list)
% rho2Mは上のfragmentに関して、第2qubitをM（X, Y, Z）
% の基底(|0>, |+>など)で入力した後の第2,3qubitの密度演算子の出力
% Vは最初のfragmentの演算子（演算子の構成に注意）

rho2M = zeros(4,4);
for j=1:2
    Mj = M_list(:,(2*j-1):2*j); 
    Mj_0 = kron(Mj,[1 0; 0 0]);
    V_Mj_0_V = V * Mj_0 * V';
    rho2M = rho2M + (3-2*j) * V_Mj_0_V;
end
end