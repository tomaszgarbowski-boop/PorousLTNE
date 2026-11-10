function v = psat(T)
%PSAT Buck saturation pressure over liquid water [Pa], T in Celsius.
v=611.21.*exp((18.678-T./234.5).*T./(257.14+T));
end
