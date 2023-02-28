function [F, G] = FMKKM_Opt_fun(X,  A, C)
%F(X) = -0.5*Tr(X'*A*X)
  G = -2*(A*X)-C;
  F = -1 *trace(X'*A*X + X'*C);
end