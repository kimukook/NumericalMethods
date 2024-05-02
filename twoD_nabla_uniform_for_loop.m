function [NABx,NABy]=twoD_nabla_uniform_for_loop(Nxp,Nyp,delta)
% Define the 2D nabla operator using 2nd-order FD methods on a uniform Nxp x Nyp grid
% Note the care done at the boundaries and corners
for i=1:Nxp
    for j=1:Nyp
        % First identify center, north, south, east, and west points
        c=i+(j-1)*Nxp; n  =c+Nxp;   s  =c-Nxp;   e  =c+1; w  =c-1;
                       nn =c+2*Nxp; ss =c-2*Nxp; ee =c+2; ww =c-2; 
        % The edge operators are defined via second-order forward/backward differences into the domain
        if     i==1;   NABx(c,c)=-3/2; NABx(c,e)= 2; NABx(c,ee)=-1/2; % West edge
        elseif i==Nxp; NABx(c,c)= 3/2; NABx(c,w)=-2; NABx(c,ww)= 1/2; % East edge
        else           NABx(c,e)= 1/2; NABx(c,w)=-1/2;                % Interior
        end
        if     j==1;   NABy(c,c)=-3/2; NABy(c,n)= 2; NABy(c,nn)=-1/2; % South edge
        elseif j==Nyp; NABy(c,c)= 3/2; NABy(c,s)=-2; NABy(c,ss)= 1/2; % North edge
        else           NABy(c,n)= 1/2; NABy(c,s)=-1/2;                % Interior
        end
    end % for j
end     % for i
NABx=NABx/delta; NABy=NABy/delta; % apply scaling
end % function NR_2D_nabla_uniform
