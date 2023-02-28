function [F, G] = FMKKM_Opt_fun2(X,  A, B, C)
%F(X) = -Tr(X'*A*X + B*X*C)
  G = -2*(A*X)-B'*C';
  F = -1 *trace(X'*A*X + B*X*C);
end