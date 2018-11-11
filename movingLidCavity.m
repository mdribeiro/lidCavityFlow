% SIMPLE Algorithm for the solution of the lid-cavity problem
% Discretization scheme used: Hybrid-scheme
%
% Mateus Dias Ribeiro
% 16.03.2014

clear all;
clc;

% Problem variables
L = 1; %input('Enter the cavity size L: [m] \n');
U = 0.1; %input('Enter the velocity U: [m/s] \n');
mi = 1; %input('Enter the viscosity: [Ns/m2] \n');
ro = 1; %input('Enter the density: [kg/m3] \n');
alpha = 0.5; % input('Enter the under-relaxation factor alpha \n');
n = 19; %input('Enter the grid order n: \n');

gamma = mi;
Re = ro*U*L/mi;

% Definitions of the u, v and P fields
u0 = zeros(n);
u = zeros(n);
v0 = zeros(n);
v = zeros(n);
P0 = zeros(n);
Pl = zeros(n);
P = zeros(n);

x = linspace(0,1,n);
y = x;

deltax = L/n;
deltay = L/n;

Ax = deltay*1;
Ay = deltax*1;

Dx = gamma*Ax/deltax;
Dy = gamma*Ay/deltay;

% Inference of the U value on the top of the cavity
for(i=1:n)
u0(1,i) = U;
u(1,i) = U;
end

res = 10;
st=0;
cont=0;

while st~=1
    % Obtaining u0 field
in =1;
while in<=30

for(i=2:n-1)
    for(j=2:n-1)
        Fe = 0.5*ro*Ax*(u(i,j)+u(i,j+1));
        Fw = 0.5*ro*Ax*(u(i,j)+u(i,j-1));
        Fn = 0.5*ro*Ay*(v(i-1,j)+v(i-1,j-1));
        Fs = 0.5*ro*Ay*(v(i,j)+v(i,j-1));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        ap = ae + aw + an + as + (Fe-Fw) + (Fn-Fs);
      
        u0(i,j) = (ae*u0(i,j+1) + aw*u0(i,j-1) + an*u0(i-1,j) + as*u0(i+1,j) + (P0(i,j-1) - P0(i,j))*Ax)/ap;
    end
end

in=in+1;
end

    % Obtaining v0 field
in =1;
while in<=30

for(i=2:n-1)
    for(j=2:n-1)
        Fe = 0.5*ro*Ax*(u(i,j+1)+u(i+1,j+1));
        Fw = 0.5*ro*Ax*(u(i,j)+u(i+1,j));
        Fn = 0.5*ro*Ay*(v(i,j)+v(i-1,j));
        Fs = 0.5*ro*Ay*(v(i,j)+v(i+1,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        ap = ae + aw + an + as + (Fe-Fw) + (Fn-Fs);
      
        v0(i,j) = (ae*v0(i,j+1) + aw*v0(i,j-1) + an*v0(i-1,j) + as*v0(i+1,j) + (P0(i+1,j) - P0(i,j))*Ay)/ap;
    end
end

in=in+1;
end

% Obtaining Pl (Pressure deviation)
in =1;
while in<=30

for(i=1:n)
    for(j=1:n)
    if((i==1)&&(j==1))   % corner up-left
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = 0;
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=0;
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - 0)*Ay;
      
        %Pl(i,j) = (aE*Pl(i,j+1) + aS*Pl(i+1,j) + b(i,j))/aP;
        Pl(i,j) = 0;
        
    elseif((i==1)&&(j==n))   % corner up-right
        Fe = 0;
        Fw = ro*Ax*(u(i,j));
        Fn = 0;
        Fs = ro*Ay*(v(i,j));
      
        ae=0;
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=0;
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        dw = Ax/aw;
        ds = Ay/as;
      
        aW = ro*dw*Ax;
        aS = ro*ds*Ay;
      
        aP = aW + aS;
      
        b(i,j) = (ro*u0(i,j) - 0)*Ax + (ro*v0(i,j) - 0)*Ay;
      
        %Pl(i,j) = (aW*Pl(i,j-1) + aS*Pl(i+1,j) + b(i,j))/aP;
        Pl(i,j) = 0;
        
     elseif((i==n)&&(j==1))   % corner down-left
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aE*Pl(i,j+1) + aN*Pl(i-1,j) + b(i,j))/aP;
        %Pl(i,j) = 0;
        
     elseif((i==n)&&(j==n))   % corner down-right
    Fe = 0;
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=0;
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - 0)*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aW*Pl(i,j-1) + aN*Pl(i-1,j) + b(i,j))/aP;
        %Pl(i,j) = 0;
    
    elseif((i~=1)&&(i~=n)&&(j==1))   % wall left
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aE*Pl(i,j+1) + aN*Pl(i-1,j) + aS*Pl(i+1,j) + b(i,j))/aP;
        
     elseif((i~=1)&&(i~=n)&&(j==n))   % wall right
        Fe = 0;
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=0;
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - 0)*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aW*Pl(i,j-1) + aN*Pl(i-1,j) + aS*Pl(i+1,j) + b(i,j))/aP;
        
     elseif((j~=1)&&(j~=n)&&(i==1))   % wall top
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = 0;
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=0;
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - 0)*Ay;
      
        %Pl(i,j) = (aE*Pl(i,j+1) + aW*Pl(i,j-1) + aS*Pl(i+1,j) + b(i,j))/aP;
        Pl(i,j) = 0;
        
     elseif((j~=1)&&(j~=n)&&(i==n))   % wall bottom
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aE*Pl(i,j+1) + aW*Pl(i,j-1) + aN*Pl(i-1,j) + b(i,j))/aP;
        %Pl(i,j) = 0;
     
     else     
        Fe = ro*Ax*(u(i,j+1));
        Fw = ro*Ax*(u(i,j));
        Fn = ro*Ay*(v(i-1,j));
        Fs = ro*Ay*(v(i,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
      
        de = Ax/ae;
        dw = Ax/aw;
        dn = Ay/an;
        ds = Ay/as;
      
        aE = ro*de*Ax;
        aW = ro*dw*Ax;
        aN = ro*dn*Ay;
        aS = ro*ds*Ay;
      
        aP = aE + aW + aN + aS;
      
        b(i,j) = (ro*u0(i,j) - ro*u0(i,j+1))*Ax + (ro*v0(i,j) - ro*v0(i-1,j))*Ay;
      
        Pl(i,j) = (aE*Pl(i,j+1) + aW*Pl(i,j-1) + aN*Pl(i-1,j) + aS*Pl(i+1,j) + b(i,j))/aP;
     end
    end
end

in=in+1;
end

% P correction
P = P0 + alpha*Pl;
P0 = P;


% u, v corrections
for(i=2:n-1)
    for(j=2:n-1)
        Fe = 0.5*ro*Ax*(u(i,j)+u(i,j+1));
        Fw = 0.5*ro*Ax*(u(i,j)+u(i,j-1));
        Fn = 0.5*ro*Ay*(v(i-1,j)+v(i-1,j-1));
        Fs = 0.5*ro*Ay*(v(i,j)+v(i,j-1));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
        ap = ae + aw + an + as;
      
        dx = Ax/ap;
            
        u(i,j) = u0(i,j) + dx*(Pl(i,j-1) - Pl(i,j));
        u(i,j) = alpha*u(i,j) + (1-alpha)*u0(i,j);
    end
end
for(i=2:n-1)
    for(j=2:n-1)
        Fe = 0.5*ro*Ax*(u(i,j+1)+u(i+1,j+1));
        Fw = 0.5*ro*Ax*(u(i,j)+u(i+1,j));
        Fn = 0.5*ro*Ay*(v(i,j)+v(i-1,j));
        Fs = 0.5*ro*Ay*(v(i,j)+v(i+1,j));
      
        ae=abs(max(max(-Fe,(Dx-Fe/2)),0));
        aw=abs(max(max(Fw,(Dx+Fw/2)),0));
        an=abs(max(max(-Fn,(Dy-Fn/2)),0));
        as=abs(max(max(Fs,(Dy+Fs/2)),0));
        ap = ae + aw + an + as;
      
        dy = Ay/ap;
            
        v(i,j) = v0(i,j) + dy*(Pl(i+1,j) - Pl(i,j));
        v(i,j) = alpha*v(i,j) + (1-alpha)*v0(i,j);
    end
end


res=(max(sqrt(mean((Pl).^2))))
res2=(max(sqrt(mean((u-u0).^2))))
res3=(max(sqrt(mean((v-v0).^2))))
rms1(cont+1) = res;
rms2(cont+1) = res2;
rms3(cont+1) = res3;

if res<0.0001
%if cont==1000
    st=1;
end

cont = cont + 1;

u0 = u;
v0 = v;

contourf(x,y,u,20), set(gca,'Ydir','reverse');
%pause
%contourf(x,y,P,20), set(gca,'Ydir','reverse');
%pause
%contourf(x,y,U,20), set(gca,'Ydir','reverse');
%pause
%contourf(x,y,P,20), set(gca,'Ydir','reverse');
%h = streamslice(x,y,u,-v); set(gca,'Ydir','reverse');
%set(h,'color','black');
drawnow;

end

pause


% Finding the center of the vortex
N = n*n;
%ip = 5;
%jp = 11;
for(ip=1:n)
    for(jp=1:n)

for (i=1:n)
    for (j=1:n)
        dx = j - jp;
        dy = ip - i;
        U = [u(i,j) v(i,j) 0];
        PM = [dx dy 0];
        angle(i,j) = atan2(norm(cross(U,PM)),dot(U,PM));
        sinangle(i,j) = sin(angle(i,j));
    end
end

I(ip,jp) = sum(sum(sinangle));
gamma1(ip,jp) = (1/N)*I(ip,jp);

    end

end

contourf(x,y,gamma1,20), set(gca,'Ydir','reverse');

% Calculating U Magnitude
for (i=1:n)
    for (j=1:n)
        U(i,j) = sqrt(u(i,j)^2+v(i,j)^2);
    end
end

pause
contourf(x,y,P,20), set(gca,'Ydir','reverse');
pause
contourf(x,y,U,20), set(gca,'Ydir','reverse');


        
