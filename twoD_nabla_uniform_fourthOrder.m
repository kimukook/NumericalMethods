function [NABx,NABy]=twoD_nabla_uniform_fourthOrder(Nx,Ny,delta)
% Define the 2D nabla operator using 4th-order FD methods on a uniform Nx x Ny grid
% Note the care exercised at the boundaries.
for i=1:Nx
 for j=1:Ny
 % First identify center, north, south, east, and west points
 c=i+(j-1)*Nx; n =c+Nx; s =c-Nx; e =c+1; w =c-1;
 nn =c+2*Nx; ss =c-2*Nx; ee =c+2; ww =c-2; 
 nnn =c+3*Nx; sss =c-3*Nx; eee =c+3; www =c-3;
 nnnn=c+4*Nx; ssss=c-4*Nx; eeee=c+4; wwww=c-4;
 if i==1; NABx(c,c)=-25/12; NABx(c,e)= 4; NABx(c,ee)=-3; NABx(c,eee)= 4/3; NABx(c,eeee)=-1/4; % West edge
 elseif i==Nx; NABx(c,c)= 25/12; NABx(c,w)=-4; NABx(c,ww)= 3; NABx(c,www)=-4/3; NABx(c,wwww)= 1/4; % East edge
 elseif i==2; NABx(c,w)=-1/4; NABx(c,c)=-5/6; NABx(c,e) = 3/2; NABx(c,ee) =-1/2; NABx(c,eee) = 1/12;
 elseif i==Nx-1; NABx(c,e)= 1/4; NABx(c,c)= 5/6; NABx(c,w) =-3/2; NABx(c,ww) = 1/2; NABx(c,www) =-1/12;
 else NABx(c,ww)=1/12; NABx(c,w)=-8/12; NABx(c,e)=8/12; NABx(c,ee)=-1/12; % Interior
 end % if i
 if j==1; NABy(c,c)=-25/12; NABy(c,n)= 4; NABy(c,nn)=-3; NABy(c,nnn)= 4/3; NABy(c,nnnn)=-1/4; % South edge
 elseif j==Ny; NABy(c,c)= 25/12; NABy(c,s)=-4; NABy(c,ss)= 3; NABy(c,sss)=-4/3; NABy(c,ssss)= 1/4; % North edge
 elseif j==2; NABy(c,s)=-1/4; NABy(c,c)=-5/6; NABy(c,n) = 3/2; NABy(c,nn) =-1/2; NABy(c,nnn) = 1/12;
 elseif j==Ny-1; NABy(c,n)= 1/4; NABy(c,c)= 5/6; NABy(c,s) =-3/2; NABy(c,ss) = 1/2; NABy(c,sss) =-1/12;
 else NABy(c,ss)=1/12; NABy(c,s)=-8/12; NABy(c,n) = 8/12; NABy(c,nn) =-1/12; % Interior
 end % if j
 end % for j
end % for i
NABx=NABx/delta; NABy=NABy/delta; % apply scaling
end % function NR_2D_nabla_uniform
