function B = angleNormalized(A)
%ANGLENORMALIZED Summary of this function goes here
%   Detailed explanation goes here
B = A;
B(B < 0) = B(B < 0) + 2*pi;
B(B > 2*pi) = B(B > 2*pi) - 2*pi;

end
