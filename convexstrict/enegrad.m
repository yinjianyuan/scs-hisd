function [Grad,ene] = enegrad(x,s)

d = length(x);
s = zeros(d,1); a = s;
for i=1:d
    a(i) = 5*i-4;
    if i<=5
        s(i) = -1;
    else
        s(i) = 1;
    end
end
ene = 0;
for i=1:d
    ene = ene + s(i)*a(i)*(exp(x(i))-x(i))/10;
end

Grad = (s.*a).*(exp(x)-1)/10;
Grad = -Grad;
end