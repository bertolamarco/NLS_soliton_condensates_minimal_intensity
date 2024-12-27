function Traj=QTrajectory(Q, SpecialPoints, Length,Startpoint, epsilon)
%% X=QTrajectory(Q, SpecialPoints, Length,Startpoint, epsilon,flag)
%% Computes the horizontal trajectory of the quadratic differential Q =Q(z) dz^2
%% through Startpoint, of approximate length Length, with step epsilon. 
%% The trajectory stops at SpecialPoints if encountered.
%
% August 2024
% (c) Marco Bertola
% email: Marco.Bertola@concordia.ca

L = Length/epsilon;
N=floor(L);
%N should be the number of steps so that the final length of the trajectory
%is what we requested in Length, given our adaptive step.

X=zeros(2,N);
Directions=X;
T=cell(1,2);

for leg=1:2
    X(leg,1)= Startpoint;
    x=X(leg,1);
    %y = sqrt( polyval(Num,x)/polyval(Den,x));
    y=sqrt(Q(x));
    Direction = conj(y)/abs(y);
    %ZEROPOLES= cat(2,roots(Num), roots(Den));
    %ZEROPOLES=SpecialPoints; %These are the points where trajectories should end.
    %Wrong name; anyways, it contains all zeros and poles
    %remlist= find();

    % remlist = [];
    % for jt = 1:length(ZEROS)
    %     %All critical points closer than 5(3)epsilon from the starting one are
    %     %excluded from the proximity test (because presumably have been already
    %     %counted...)
    %     if (abs(ZEROS(jt) - Startpoint)<3*epsilon)
    %        remlist=[remlist,jt];
    %     end
    % end
    %tmpZEROS=ZEROPOLES(abs(ZEROPOLES- Startpoint)>=3*epsilon).'; % We will check proximity only with these.
    
    
    if (leg==2) 
        Direction=-Direction; %We find the leg in the other direction!
    end
    Directions(leg,1)=Direction;
  
    for kt=2:N
        dt=epsilon;
    
        K1 = Direction;
        %ynew =  sqrt( polyval(Num,x+dt*K1/2)/polyval(Den,x+dt*K1/2));
        ynew = sqrt(Q(x+dt*K1/2));
        DirectionNew = conj(ynew)/abs(ynew);
        if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction)) %Continuity check on the determination of the root.
            Direction1=-DirectionNew;
        else
            Direction1=DirectionNew;
        end
    
        K2 = Direction1;
        %ynew =  sqrt( polyval(Num,x+dt*K2/2)/polyval(Den,x+dt*K2/2));
        ynew = sqrt(Q(x+dt*K2/2));
    
        DirectionNew = conj(ynew)/abs(ynew);
        if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction))
            Direction2=-DirectionNew;
        else
            Direction2=DirectionNew;
        end
    
        K3 = Direction2;
        %ynew =  sqrt( polyval(Num,x+dt*K3)/polyval(Den,x+dt*K3));
        ynew = sqrt(Q(x+dt*K3));
        DirectionNew = conj(ynew)/abs(ynew);
        if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction))
            Direction3=-DirectionNew;
        else
            Direction3=DirectionNew;
        end
        K4=Direction3;
    
        Direction= 1/6*(K1+2*K2+2*K3+K4);
        x= x+dt*Direction;
        Direction= Direction/abs(Direction);
        % Runge Kutta 4
        % End Runge Kutta

        Directions(leg,kt)=Direction;       
        X(leg,kt) = x;
    end
end
%Now we find out if our trajectory came close to some of the other points,
%in which case we truncate.

for leg=1:2
    numpoints=N;
    ProximityRadius = max(epsilon, 1e-3);
    tmp=length(SpecialPoints);
    tmp2 = abs((X(leg,:).')*ones(1,tmp)  - ones(N,1)*SpecialPoints); %Contains the distances of trajectory X from each of the zeros/poles
    %We now find the first occurrence when X gets closer than ProximityRadius
    %to any of the zeros/poles.

    [row,col]=find (tmp2<ProximityRadius); %Find indices where the trajectory is near a special point
    flags = real(Directions(leg, row).*conj(X(leg,row)-SpecialPoints(col)))<0;


    % for j=1:size(row,1)
    %     if real(Directions(leg, row(j))*conj(X(leg,row(j))-SpecialPoints(col(j))))<0 %We are pointing towards a special point! STOP!
    %         flags(j)=1;
    %     end
    % end
    if size(row(flags==1),1)>0
        numpoints = min(N,min(row(flags==1))); %Find the first instance. The tail after that we truncate below.
    end
    T{leg} = X(leg,1:numpoints);
end
if size(T{1},2)>size(T{2},2)
    Traj=[T{2}(end:-1:1), T{1}];
else
    Traj=[T{1}(end:-1:1), T{2}];
end

end
