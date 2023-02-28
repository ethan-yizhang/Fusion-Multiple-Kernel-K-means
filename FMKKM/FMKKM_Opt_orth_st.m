function [X,F]= FMKKM_Opt_orth_st(X, fun, opts, varargin)


%% Size information
if isempty(X)
    error('input X is an empty matrix');
else
    [n, k] = size(X);
end

if nargin < 2; error('[X, out]= OptStiefelGBB(X0, @fun, opts)'); end
if nargin < 3; opts = [];   end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'tau');       opts.tau  = 5; end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.2; end
if ~isfield(opts, 'retr');      opts.retr = 0; end
if ~isfield(opts, 'tiny');      opts.tiny = 1e-13; end

% copy parameters
tau     = opts.tau;
rhols   = opts.rhols;
eta     = opts.eta;
retr    = opts.retr;
tiny    = 1e-13;
%-------------------------------------------------------------------------------

%% Initial function value and gradient
% prepare 
[F,  G] = feval(fun, X , varargin{:}); 
GX = G'*X;

if retr == 1
    invH = true; if k < n/2; invH = false;  eye2k = eye(2*k); end
    if invH
        GXT = G*X';  H = 0.5*(GXT - GXT');  RX = H*X;
    else
        U =  [G, X];    V = [X, -G];       VU = V'*U;
        %U =  [G, X];    VU = [GX', X'*X; -(G'*G), -GX];
        %VX = VU(:,k+1:end); %VX = V'*X;
        VX = V'*X;
    end
end
dtX = G - X*GX;     nrmG  = norm(dtX, 'fro');
  

    XP = X;     FP = F;  
     % scale step size
    nls = 1; deriv = rhols*nrmG^2; %deriv
    while 1
        % calculate G, F,        
        if retr == 1
            if invH
                [X, ~] = linsolve(eye(n) + tau*H, XP - tau*RX);
            else
                [aa, ~] = linsolve(eye2k + (0.5*tau)*VU, VX);
                X = XP - U*(tau*aa);
            end
        else
            [X, ~] = myQR(XP - tau*dtX, k);
        end
        
        if norm(X'*X - eye(k),'fro') > tiny; X = myQR(X,k); end
        
        [F,G] = feval(fun, X, varargin{:});
        
        if F <= FP - tau*deriv || nls >= 8
            break;
        end
        tau = eta*tau;          nls = nls+1;
    end  
    %% 判断F是否下降（tau已足够小，F仍未下降，则停止迭代并撤销该步）
    if  F >= FP
        X = XP;
    end
    



end


function [Q, RR] = myQR(XX,k)
[Q, RR] = qr(XX, 0);
diagRR = sign(diag(RR)); ndr = diagRR < 0;
if nnz(ndr) > 0
    Q = Q*spdiags(diagRR,0,k,k);
    %Q(:,ndr) = Q(:,ndr)*(-1);
end
end

