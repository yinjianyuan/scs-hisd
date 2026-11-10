function hv = hv(x, v, s)

d = length(x);
hv = zeros(d,1);
s = zeros(d,1); a = s;
for i=1:d
    a(i) = 5*i-4;
    if i<=5
        s(i) = -1;
    else
        s(i) = 1;
    end
end
for i=1:d
    hv(i) = s(i)*a(i)*exp(x(i))/10*v(i);
end
end