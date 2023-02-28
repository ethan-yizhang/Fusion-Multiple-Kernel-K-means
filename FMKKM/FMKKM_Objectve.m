function [obj,obj1,obj2,obj3] = FMKKM_Objectve(KH,Hstar,HP,WP,Alpha,Beta,Gamma,lambda1,lambda2)
%F(X) = -0.5*Tr(X'*A*X)
  numker = size(HP,3);
  Kmatrix = sumKbeta(KH,Gamma.^2);
  obj1 = - trace(Hstar' * Kmatrix * Hstar) + trace(Kmatrix);
  obj2 = 0;
  obj3 = 0;
  for p = 1:numker
     obj2 = obj2 + Alpha(p).^2 * ( - trace(HP(:,:,p)'*KH(:,:,p)*HP(:,:,p)) + trace(KH(:,:,p)));
     obj3 = obj3 + trace(Hstar' * Beta(p) * HP(:,:,p) * WP(:,:,p));
  end
  
  obj = obj1 + lambda1 * obj2 - lambda2 * obj3;
end