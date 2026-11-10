function hv = hv(x, v, s)

d = length(x);
x = x(:);
v = v(:);

hv = zeros(d,1);

ss = ones(d,1);
ss(1:5) = s;
ss(d) = 150;

atan_term = ss .* (2 + 4*(1-x).*atan(x-1)) ./ (1+(x-1).^2).^2;

diagH = zeros(d,1);
diagH(1) = 1200*x(1)^2 - 400*x(2) + 2 + atan_term(1);

for i = 2:d-1
    diagH(i) = 1200*x(i)^2 - 400*x(i+1) + 202 + atan_term(i);
end

diagH(d) = 200 + atan_term(d);

hv = diagH .* v;

for i = 1:d-1
    c = -400*x(i);
    hv(i)   = hv(i)   + c * v(i+1);
    hv(i+1) = hv(i+1) + c * v(i);
end

end
