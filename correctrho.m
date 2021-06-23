function [correct_rho, correct_vector] = correctrho(U1,V)
% |000><000|を入力した際の理論上の最終出力としての密度演算子と状態ベクトルの算出
% U1は上のfragmentの2qubit演算, Vは後ろのfragmentの2qubit演算
% 演算子の順番に注意（U1の構成に関して、演算の順番はCX,HIなのかHI,CXなのか、等...）

% initial_stateは初期状態|000><000|(密度演算子)
initial_state = zeros(8,8);
initial_state(1,1) = 1;

I = eye(2);

% all_operatorは回路全体としての演算子
all_operator = kron(I,V) * kron(U1, I);

% 最終出力の理論上の密度演算子
correct_rho = all_operator * initial_state * all_operator';

% 最終出力の理論上の状態ベクトル
correct_vector = all_operator * [1;0;0;0;0;0;0;0];

end


