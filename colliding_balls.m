
clc 
close all
clc
clear all
%-----------------------
N=1+randi(7);                  % number of particles
L=3;  
Lx=3;
Ly=2;%length of the box
dt=0.00005; 

X=zeros(2,N); 
V=zeros(2,N); 
R = 0.1 + 0.2 * rand(1, N);
M=R.^2; % mass of each particle
theta=0:2*pi/20:2*pi;
nstep=10000*2;               % number of time steps
pit=80;                      % number of iterations between two plots
collis=0;
%initilizing random position and velocity 
for i = 1:N
    ok = false;
    while ~ok
        % Generate random position
        X(:, i) = [0.1 + 2.5 * rand(), 0.1 + 1.8 * rand()]';

        % Check if the ball is within the specified rectangle
        if X(1, i) - R(i) >= 0.3 && X(1, i) + R(i) <= 2.7 && ...
           X(2, i) - R(i) >= 0.3 && X(2, i) + R(i) <= 1.7
            % Check for collision with other balls
            flag = false;
            for j = 1:i-1
                if norm(X(:, j) - X(:, i)) <= R(i) + R(j) + 0.01
                    flag = true;
                    break;
                end
            end

            if ~flag
                ok = true; % Position is valid
            end
        end
    end
end


  Ax=zeros(N,N);
  Ay=zeros(N,N);
  Ar=zeros(N,N);
  for i=1:N
    for j=i:N
        Ar(i,j)=R(1,i)+R(1,j);
    end
  end

 hold on
 axis equal
 for i=1:N       
  x=X(1,i)+R(1,i)*cos(theta);
  y=X(2,i)+R(1,i)*sin(theta);
 
 end

 pause(1)
 randomArray = randi([0, 1], 1, N) * 2 - 1;  
 randomArray2 = 10 + 10 * rand(1, N);
 V(1,:)= randomArray.*randomArray2;
 randomArray = randi([0, 1], 1, N) * 2 - 1;  
 randomArray2 = 10 + 10 * rand(1, N);
 V(2,:)=randomArray.*randomArray2;
 



for T=0:nstep
    X=X+dt*V;
    %check if particules collided with eachother
    for i=1:N
        Ax(i,:)=X(1,i)-X(1,:);
        Ay(i,:)=X(2,i)-X(2,:);
        % I used this to Calculate kinetic energy
    end
    %Ax=triu(Ax);
    %Ay=triu(Ay);
    Nrm=(Ax.^2+Ay.^2).^(0.5)-Ar-10^-3; % calculate the distance between each two particles
    Nrm=triu(Nrm(1:N-1,2:N));
    [row,col]=find(Nrm<0); % find the particles that collided 
    l=length(row);
    %if particles collided calculate the new velocities
    if l~=0
        col=col+1;       
        for t=1:l
            i=row(t,1);
            j=col(t,1);
            C1=(2*M(1,j)/(M(1,i)+M(1,j)))*dot(V(:,i)-V(:,j),X(:,i)-X(:,j))/(norm(X(:,i)-X(:,j))^2);
            C2=(2*M(1,i)/(M(1,i)+M(1,j)))*dot(V(:,j)-V(:,i),X(:,j)-X(:,i))/(norm(X(:,i)-X(:,j))^2);
            V(:,i)=V(:,i)-C1*(X(:,i)-X(:,j));
            V(:,j)=V(:,j)-C2*(X(:,j)-X(:,i));
            collis=collis+1;

        end
    end
 
     for i=1:N
        %check if particules collided with the walls horizontally
        if X(1,i)+R(1,i)>=3-10^-3 || X(1,i)-R(1,i)<=+10^-3
            V(1,i)=-V(1,i);
            collis=collis+1;

        end
        %check if particules collided with the walls vertically
        if X(2,i)+R(1,i)>=2-10^-3 || X(2,i)-R(1,i)<=0+10^-3
            V(2,i)=-V(2,i);
            collis=collis+1;

        end
       

     end






     %plot
     if T==pit*ceil(T/pit)
         clf
         hold on
         axis equal
         for i=1:N      
              x=X(1,i)+R(1,i)*cos(theta);
              y=X(2,i)+R(1,i)*sin(theta);
              color=["b","g","k","y","c","m","b","g","k","y","c","m","b","g","k","y","c","m","b","g","k","y","c","m","b","g","k","y","c","m","b","g","k","y","c","m"];
              xlim([-1,4])
              ylim([-1,3])

              plot([0 3],[0 0],"r",[0 0],[0 2],"r",[3 3],[0 2],"r",[0 3],[2,2],"r","LineWidth",2)
              grid on
              plot(x,y,color(i))
              fill(x,y,color(i))


         end
        title(T*dt)
        text(0, 0, sprintf('Collision: %d', collis), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');        drawnow
        pause(1e-10)
     end
     
end


