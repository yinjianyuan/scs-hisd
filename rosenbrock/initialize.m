
s1 = s1s(casenumber);
d = 1000; 

x = ones(d,1);
% rng('shuffle');
% rng(1)
% rao = randn(d,1);
% x = x + 1*rao/norm(rao)

A = zeros(d,d); 

for i = 1:d
    if i <=5
        si = s1;
    elseif i == d
        si = 150;
    else
        si = 1;
    end

    atan_term = si * (2 + 4*(1-x(i))*atan(x(i)-1)) / (1+(x(i)-1)^2)^2;

    if i == 1
        A(i,i) = 1200*x(i)^2 - 400*x(i+1) + 2 + atan_term;
    elseif i == d
        A(i,i) = 200 + atan_term;
    else
        A(i,i) = 1200*x(i)^2 - 400*x(i+1) + 202 + atan_term;
    end
end

for i = 1:d-1
    A(i,i+1) = -400*x(i);
    A(i+1,i) = -400*x(i);
end

[eigenvalue]=eig(A);
disp(eigenvalue(4))
