function [x, output] = hiosd_lobpcg(F,x,options)
%% reading parameters
n = length(x); % n is the number of variables
if isfield(options,'innpfunc')
    innp = @(x,y)options.innpfunc(x,y);
    nor = @(x)sqrt(innp(x,x));
else
    disp('HiOSD: using the standard inner product.');
    innp = @(x,y)(x'*y);
    nor = @(x)sqrt(x'*x);
end

if isfield(options,'k')
    if isfield(options,'V')
        if options.k ~= size(options.V, 2)
            options.k = size(options.V, 2);
            disp(['HiOSD: k is reset as ' num2str(options.k) ', the column - number of V.']);
        end
    end
elseif isfield(options,'V')
    options.k = size(options.V, 2);
else
    options.k = 1;
    disp('HiOSD: k is reset as 1.');
end
k = options.k;

if isfield(options,'V')
    if rank(options.V) < k
        disp('HiOSD: V is random vectors.');
        V = randn(n,k);
    else
        V = options.V;
    end
else
    disp('HiOSD: V is random vectors.');
    V = randn(n,k);
end

for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end

if isfield(options,'dt')
    dt = options.dt;
else
    disp('HiOSD: dt = 0.01.');
    dt = 0.01;
end

if isfield(options,'l')
    l = options.l;
else
    l = 1e-6;
end

if isfield(options,'minl')
    minl = options.minl;
else
    minl = 1e-6;
end

if isfield(options,'epsf')
    epsf = options.epsf;
else
    disp('HiOSD: epsf = 1e-8.');
    epsf = 1e-8;
end

if isfield(options,'maxiter')
    maxiter = options.maxiter;
else
    maxiter = 1e4;
end

if isfield(options,'betat')
    betat = options.betat;
else
    betat = 1;
end
if isfield(options,'betau')
    betau = options.betau;
else
    betau = 1;
end

if ~isfield(options,'outputX')
    options.outputX = 0;
end

if options.outputX ~= 0
    X = zeros(k+2+n,ceil(maxiter/options.outputX));
    storing = 'X(:,ceil(iter/options.outputX))=[nor(f);alpha;beta;x];';
else
    storing = '';
end

if ~isfield(options,'outputp')
    options.outputp = 1000;
end

if ~isfield(options,'outputd')
    options.outputd = 0;
end

if options.outputd ~= 0 && isfield(options,'draws')
    drawing = options.draws;
else
    drawing = '';
end

if isfield(options,'tau')
    tau = options.tau;
else
    tau = 0.1;
end

if isfield(options,'precd')
    P = options.precd;
else
    P = @(x)x;
end

if isfield(options,'adding')
    adding = options.adding;
else
    adding = '';
end

%% setting environment
alpha (1 : k, 1) = 0;
H = @ (xx,vv,ll) - ( F (xx + ll * vv) - F (xx - ll * vv) ) / (2 * ll);
U = zeros(n,k);
W = zeros(n,k);
Y = zeros(n,k);
%% bb iterations
iter = 1;
Vp = zeros(n,k);
Up = zeros(n,k);
f = F (x);
norf = nor(f);
gp = zeros(n,1);
Dx = f - 2 * V * innp(V,f); Dx = Dx * dt;
while iter <= maxiter
    emp = true( 3 * k, 1 );
    g = f - 2 * V * innp(V,f);
    Dg = g - gp;
    gp = g;
    
    beta = abs ( innp(Dx,Dg) / innp(Dg,Dg) );
    beta = min ( beta, betat * dt);
    beta = max ( beta, betau * dt);
    if norf*beta > tau
        beta=tau/norf;
    end
    Dx = g * beta;
    x = x + Dx;
    for i = 1 : k
        U(:,i) = H (x, V(:,i), l);
        alpha(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alpha(i) * V(:,i) ) ;
        W(:,i) = P(W(:,i));
        W(:,i) = W(:,i) - V * innp( V, W(:,i) ) ;
        W(:,i) = W(:,i) - W(:,1:i-1) * innp( W(:,1:i-1), W(:,i) ) ;
        nrmW = nor(W(:,i)) ;
        if nrmW > 1e-8
            W(:,i) = W(:,i) / nrmW ;
            Y(:,i) = H (x, W(:,i), l);
        else
            %           Y(:,i) = H (x, W(:,i), l);
            W(:,i) = 0;
            Y(:,i) = 0;
            emp(k+i) = false;
        end
    end
    for i = 1 : k
        zeta = innp( [V W Vp(:,1:i-1)], Vp(:,i) );
        Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * zeta;
        %         Up(:,i) = Up(:,i) - [U Y Up(:,1:i-1)] * zeta;
        nrmVpi = nor(Vp(:,i));
        if nrmVpi > 1e-8
            Vp(:,i) = Vp(:,i) / nrmVpi ;
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp( [V W Vp(:,1:i-1)], Vp(:,i) );
            Vp(:,i) = Vp(:,i) / nor(Vp(:,i)) ;
            %             Up(:,i) = Up(:,i) / nrmVpi ;
            Up(:,i) = H (x, Vp(:,i), l);
        else
            %           Up(:,i) = H (x, Vp(:,i), l);
            Vp(:,i) = 0;
            Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end
    UU = [V,W,Vp]; UU = UU (:,emp);
    YY = [U,Y,Up]; YY = YY (:,emp);
    
    Pn = innp(UU,YY);  Pn = (Pn+Pn') / 2;
    [eta, alpha] = eig( Pn, 'vector' );
    
    Vp=V;
    V = UU * real(eta(:,1:k));
    alpha = real(alpha(1:k));
    
    l = max(l/(1+dt),minl);
    
    f = F (x);
    norf = nor(f);
    eval(storing);
    if mod(iter, options.outputp)==0
        disp([num2str(iter,'%d') ': nf=' num2str(norf,'%.4e') '; a=' num2str(alpha', '%.4f ')]);
    end
    if mod(iter, options.outputd)==0
        eval(drawing);drawnow;
    end
    if norf < epsf && dt>0
        if options.outputX~=0
            X = X(:,1:ceil(iter/options.outputX));
        end
        break;
    end
    if mod(iter,10)==0
        grad(round(iter/10)) = norf;
    end
    iter = iter+1;
    eval(adding);
end
if iter==maxiter+1
    iter=maxiter;
end

if dt == 0
    x = V;
    output = alpha;
    return
end

output = struct('x', x, 'V', V, 'it', iter, 'al', alpha,'gr',grad);
if options.outputX ~= 0
    output.X = X;
end

end