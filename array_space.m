function [a_s] = array_space(M,N,d,lambda,theta,phi)
    a_z = exp(1j*2*pi/lambda*(0:M-1)*d*sind(phi)).';
    a_y = exp(-1j*2*pi/lambda*(0:N-1)*d*cosd(phi)*sind(theta)).';
    a_s = kron(a_y,a_z);
end