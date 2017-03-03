function [p] = convAdvec(h,sigma,mu,q)

%   Input:
%           h: scalar. cellsize
%       sigma: vector. spatially varying coefficient, evaluated at cell centres.
%          mu: scalar. constant coefficient
%           q: vector. Source function, evaluated at nodes.
%  
%   Output: 
%           p: Approximate value of p on the nodes.
%           
%   You must write your discretization as a matrix equation A*p = q.
%   Once you have formed the matrix A, you can solve for p using the 
%   Matlab backslash operator. Make sure that you form A as a sparse
%   matrix. If A is dense the backslash operator will be very inefficient.
 %% Creating centers, nodes and the inside centers
n=1/h;
xN = [0:h:1];
xI = xN(2:n);
xC = [(h/2):h:(1-(h/2))];
%% Creating Diffusion Portion
Dnc=spdiags(ones((n),1)*[-1 1],[0,1],(n)-1,(n))/h;
Dcn = -transpose(Dnc);
diff = Dnc*sigma*Dcn;
%% Creating Convection
Avg = spdiags(ones((n),1)*[1 1],[0,1],(n)-1,(n))/2;
conv = mu.*Avg*Dcn;
%% Creating A
A = diff + conv;
%% Creating P
p = A\q'


% %%Test Function of ConvAdvec to produce sin graph. experimental and
% %%theoretical both plotted at the end
% %% used different file to test the code above and it worked :D 

% for i = 2:4
% h=10.^(-i);
% n=1/h;
% % x on the nodes
% xN = [0:h:1];
% % x on the nodes excluding the end points
% xI = xN(2:n);
% % x on the cell centers
% xC = [(h/2):h:(1-(h/2))];
% mu = 0.1;
% %mu = 10;
% % creating sigma points
% sig = (1+xC.^2)';
% Sig = spdiags(sig, [0], n, n);
% % q flow 
% diff = 4*pi*xI.*cos(2*pi*xI)-4*pi^2*(xI.^2 + 1).*sin(2*pi*xI);
% flow = mu*2*pi*cos(2*pi*xI);
% q =  diff + flow;
% Pexp = convAdvec(h,Sig,mu,q);
% Ptheo = sin(2*pi*xC);
% i+1;
% 
% figure
% hold on
% plot(Pexp)
% plot(Ptheo)
% hold off
% pause(10);
% end
% % %%this code stop running at h = 10^-5 and errors itself (ie when h is less
% % %%%than 10^-4
end 
