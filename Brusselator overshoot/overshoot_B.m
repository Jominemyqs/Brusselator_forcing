function b = overshoot_B(t, b0, bmax, Tup, Thold, Tdown)
%OVERSHOOT_B Piecewise linear overshoot forcing for Brusselator parameter b.
%
% b(t):
%   ramp from b0 to bmax over Tup
%   hold at bmax for Thold
%   ramp back to b0 over Tdown
%   stay at b0 afterwards

T1 = Tup;
T2 = Tup + Thold;
T3 = Tup + Thold + Tdown;

if t < 0
    b = b0;
elseif t < T1
    b = b0 + (bmax - b0) * (t / Tup);
elseif t < T2
    b = bmax;
elseif t < T3
    b = bmax - (bmax - b0) * ((t - T2) / Tdown);
else
    b = b0;
end

end