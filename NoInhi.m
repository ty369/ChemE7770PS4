function f=NoInhi(B)
% function for 7770 ps4 #2
% B(1)=[A]; B(2)=[B]; B(3)=[C];
f=[(5.*B(1))./(5+B(1))-B(2)./(5+B(2));
(5.*B(1))./(5+B(1))-B(3)./(5+B(3));
100-B(1)-B(2)-B(3)];
end
