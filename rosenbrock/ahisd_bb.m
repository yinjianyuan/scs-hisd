function [x, output]=ahisd_bb(F, H, x, options)
d = 1000; s = options.s;

iter0=1;
n = d;
if isfield(options,'innpfunc')
    innp = @(x,y)options.innpfunc(x,y);
    nor = @(x)sqrt(innp(x,x));
else
    disp('HiOSD: standard inner product.');
    innp = @(x,y)(x'*y);
    nor = @(x)sqrt(x'*x);
end

if isfield(options,'precd')
    P= options.precd;
else
    P=@(x)x;
end

if isfield(options,'gamma')
    gamma = options.gamma;
elseif isfield(options,'hb_gamma')
    gamma = options.hb_gamma;
else
    gamma = 0.5;
    disp('HiSD-HB: gamma = 0.5.');
end
if gamma < 0 || gamma >= 1
    error('HiSD-HB: options.gamma should satisfy 0 <= gamma < 1.');
end
disp(['gamma=' num2str(gamma)]);
rao = options.rao;
options.k = size(options.V,2);
k = options.k;
V = options.V;
for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end

dt=options.dt;
epsf=options.epsf;
maxiter=options.maxiter;

%% setting environment
alph=zeros(k, 1) ;
U = zeros(n,k);
W = zeros(n,k);
Y = zeros(n,k);
%% bb iterations
Vp = zeros(n,k);
Up = zeros(n,k);
[f,~] = F (x,s);
gp = zeros(n,1);
rao = options.rao;
Data(1,1)=0; Data(1,2)=rao;
xold = x; % x_{-1}=x_0 for the first heavy-ball step
Dx = f - 2 * V * innp(V,f); Dx = Dx * dt;
for iter=1:maxiter
    emp = true( 3 * k, 1 );
    [f,~] = F(x,s);

    tmpp = 2 * V * innp(V,f);
    g = f - tmpp;
    Dg = g - gp;
    gp = g;
    
    beta = abs ( innp(Dx,Dg) / innp(Dg,Dg) );
    beta = min ( beta, options.betat * dt);
    beta = max ( beta, options.betau * dt); 
    Dx = g * beta + gamma * (x - xold);
    xold = x;
    x = x + Dx;
    
    
    for i = 1 : k
        U(:,i) = H (x, V(:,i), s);
        alph(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alph(i) * V(:,i) ) ;
        W(:,i) = P(W(:,i)); 
        W(:,i) = W(:,i) - V * innp( V, W(:,i) ) ;
        if isfield(options,'orthV') 
            W(:,i) = W(:,i) - options.orthV * innp( options.orthV, W(:,i) ) ;
        end
        W(:,i) = W(:,i) - W(:,1:i-1) * innp( W(:,1:i-1), W(:,i) ) ;
        
        nrmW = nor(W(:,i)) ;
        if nrmW > 1e-8
            W(:,i) = W(:,i) / nrmW ;
            Y(:,i) = H (x, W(:,i),s);
        else
            W(:,i) = 0;
            Y(:,i) = 0;
            emp(k+i) = false;
        end
    end
    for i = 1 : k
        zeta = innp( [V W Vp(:,1:i-1)], Vp(:,i) );
        Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * zeta;
        if isfield(options,'orthV') 
            Vp(:,i) = Vp(:,i) - options.orthV * innp( options.orthV, Vp(:,i) ) ;
        end
        nrmVpi = nor(Vp(:,i));
        if nrmVpi > 1e-8
            Vp(:,i) = Vp(:,i) / nrmVpi ;
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp( [V W Vp(:,1:i-1)], Vp(:,i) );
            Vp(:,i) = Vp(:,i) / nor(Vp(:,i)) ;
            Up(:,i) = H (x, Vp(:,i),s);
        else
            Vp(:,i) = 0;
            Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end
    UU = [V,W,Vp]; UU = UU (:,emp);
    YY = [U,Y,Up]; YY = YY (:,emp);
    
    Pn = innp(UU,YY);  Pn = (Pn+Pn') / 2;
    [eta, alph] = eig( Pn, 'vector' );
    
    Vp=V;
    V = UU * real(eta(:,1:k));
    alph = real(alph(1:k));
    
    [f,~] = F (x,s);
    res = norm(x-ones(d,1));
    for i = 1 : k
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
        V(:,i) = V(:,i)/nor(V(:,i));
    end
        
    if res < epsf
        break;
    end
    if mod(iter,options.outputp)==0
        [g1,e1] = F(x,s);
        disp([num2str(iter) '   ' num2str(res) '   ' num2str(alph') '   ' num2str(norm(g1)/sqrt(n)) '   ' num2str(e1) ]);
    end
    if mod(iter,10) == 0
        iter0 = round(iter/10)+1;
        Data(iter0,1) = iter;
        Data(iter0,2) = norm(x-ones(d,1));
    end

    iter0 = iter0+1;
end
    if mod(iter,options.outputd)==0
        eval(options.draws);
    end
 
output = struct('x', x, 'V', V, 'it', iter, 'alph', alph, 'data', Data);

end
