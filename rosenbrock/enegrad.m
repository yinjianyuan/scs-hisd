function [Grad, ene] = enegrad(x, s)

% Energy and negative gradient for
% E(x) = sum_{i=1}^{d-1} [100*(x_{i+1}-x_i^2)^2 + (1-x_i)^2]
%      + sum_{i=1}^{d} s_i*atan(x_i-1)^2,
%
% Output Grad = -grad E, consistent with HiSD code using F = negative gradient.

x = x(:);
d = length(x);

% construct s_i
ss = ones(d,1);
ss(1:5) = s;
ss(d) = 150;

%% energy
ene = 0;

% Rosenbrock part
for i = 1:d-1
    ene = ene + 100*(x(i+1) - x(i)^2)^2 + (1 - x(i))^2;
end

% modified arctan part
for i = 1:d
    ene = ene + ss(i)*atan(x(i)-1)^2;
end

%% gradient of energy
gradE = zeros(d,1);

% arctan derivative
z = x - 1;
grad_atan = 2*ss.*atan(z)./(1+z.^2);

if d == 1
    gradE(1) = grad_atan(1);
else
    % i = 1
    gradE(1) = -400*x(1)*(x(2)-x(1)^2) + 2*(x(1)-1) + grad_atan(1);

    % 2 <= i <= d-1
    for i = 2:d-1
        gradE(i) = ...
            200*(x(i) - x(i-1)^2) ...
            - 400*x(i)*(x(i+1)-x(i)^2) ...
            + 2*(x(i)-1) ...
            + grad_atan(i);
    end

    % i = d
    % No term 2*(x(d)-1), because (1-x_i)^2 only goes to d-1.
    gradE(d) = 200*(x(d) - x(d-1)^2) + grad_atan(d);
end

Grad = -gradE;

end