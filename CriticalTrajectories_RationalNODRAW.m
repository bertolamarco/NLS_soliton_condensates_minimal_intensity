function [CritTraj,Adjacency,TypeCritpts, Q] = CriticalTrajectories_RationalNODRAW(varargin)
%% CriticalTrajectories_Rational(Const, Zeros, mu, Poles,nu,epsilon,L, varargin)
%
%% August 2024
%% (c) Marco Bertola
%% email: Marco.Bertola@concordia.ca
%% Usage 
%% Crittraj (Const, Z, mu, Poles, nu, epsilon,Length,interactive,Points)
%% Plots the critical trajectories of the quadratic differential
%% R(x) dx^2 where R(x) is a rational function.
%% Usage Crittraj(Const, Zeros, mu, Poles, nu);
%% Zeros is a vector of complex numbers containing the position of the zeroes;\n\
%% mu is a vector of integers containing the multiplicities of the zeroes
%%  Poles, nu are  similar (can be empty).
%% Const is the overall multiplicative constant that multiplies 
%% the rational function with those zeroes and poles.
%% epsilon is the step for the numerical integration
%% Length is the length of the arc drawn
%% interactive is a flag (0 or 1); if 1 then it prompts for points where to draw horizontal trajectories.
%% The vector Points contains points through which we want to draw
%% additional horizontal trajectories.
%% The parameter Thickness is used to scale the thickness of the lines by
%% this common factor.
%% The function returns CritTraj, a cell of critical trajectories,
%% Adjacency, the adjacency matrix of the critical points
%% TypeCritpts the order of the critical points indexed in Adjacency
%% Q the handle to the function representing the quadratic differential

p = inputParser;
CritTraj=cell(1,100); % We assume there won't be more than 100 critical trajectories! Change as needed.
countcrit=0;


addOptional(p,'Const',-1);
addOptional(p,'Zeros',[1i,-1i]);
addOptional(p,'mu',[2,2]);
addOptional(p,'Poles',[-1,1]);
addOptional(p,'nu',[1,1]);
addOptional(p,'epsilon',0.04);
addOptional(p,'L',4);
parse(p, varargin{:});
Const=p.Results.Const;
Zeros= p.Results.Zeros;
Poles= p.Results.Poles;
mu=p.Results.mu;
nu=p.Results.nu;
L=p.Results.L;
epsilon=p.Results.epsilon;

%global Den Num  Zeros_Poles ;
Length=L;


%%%%%%
%This is just really a front-end for the external function QTrajectory,
%which actually computes the horizonal trajectroy of a quadratic
%differential, but does none of the plotting or interface.

%Adjust multiplicities: if two roots are closer than epsilon to each other
%we merge them adding the corresponding multiplicities.


tmp0=length(Zeros);
tmp = ones(tmp0,1)*Zeros- Zeros.'*ones(1, tmp0);
tmp=abs(tmp) + eye(tmp0);
mindist= min(min(tmp));
while(mindist<2*epsilon)
    [t1,t2]=min(tmp);
    [~,tt2]=min(t1);
    %Zeros(t2(tt2)) and Zeros(tt2) are closer than epsilon; we need to merge
    Zeros(tt2)=[];
    mu(t2(tt2)) = mu(t2(tt2))+mu(tt2);
    mu(tt2)=[];
    tmp0=tmp0-1;
    tmp = ones(tmp0,1)*Zeros- Zeros.'*ones(1, tmp0) + eye(tmp0);
    tmp=abs(tmp);
    mindist= min(min(tmp));
end

%We should do the same with poles but it is a rare instance.. oh well
SpecialPoints = cat(2,Zeros,Poles);
ZZ= zeros(1,sum(mu));
for i=1:length(mu)
  for j=1:mu(i)
    ZZ(sum(mu(1:i-1))+j) = Zeros(i);
  end
end
PP = zeros(1,sum(nu));
for i=1:length(nu)
  for j=1:nu(i)
    PP(sum(nu(1:i-1))+j) = Poles(i);
  end
end
Num = Const * poly(ZZ);
Den = poly  (PP);
Q= @(z) polyval(Num,z)./polyval(Den,z);

%Now we compute the critical trajectories issuing from the
%simple poles   
for k=1:length(Poles) %We also plot the critical trajectories from simple poles of the quadratic differential
  if (nu(k)==1) 
    tmp=polyder(Den);
    tmp1 = (polyval(Num,Poles(k))/polyval(tmp,Poles(k)));
    tmp2 = tmp1;
    tmp2 = 1/(tmp2/abs(tmp2));
    StartAngle = tmp2;
    Phase = StartAngle;
    pt = Poles(k) + epsilon/9*Phase;
    Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
    countcrit=countcrit+1;
    CritTraj{countcrit} = Y;
  end
end

CritPoints=[Zeros, Poles(nu==1)]; % the critical points are the zeros and poles of first order.
TypeCritpts= [mu, nu(nu==1)];
Adjacency = zeros(length(CritPoints) + 1); 
%this is the adjacency matrix telling us which critical points are
%connected to which other, and by which trajectory. The last row/column is
%for "hanging'' trajectories

for k=1:length(Zeros)
    h=1e-3;
    C=Q(Zeros(k)+h)/h^(mu(k));
    StartAngle=exp(1i*angle(C)/2);

    for j=0:(mu(k)+2)-1
        Phase = StartAngle*exp(2i*pi/(mu(k)+2)*(j));
        pt = Zeros(k) + epsilon/50*Phase;
        Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
        countcrit=countcrit+1;
        CritTraj{countcrit} = Y;
        %We finally plot the trajectory in the appropriate color: thick red for
        %cuts, thick blue for (critica) shorelines, thin blue for
    end

end
%What remains to do now is to find if this trajectory connects two
%critical points...

for j=1:countcrit
    [delta_s,start] = min(abs(CritPoints - CritTraj{j}(1)));
    [delta_f, finish] = min(abs(CritPoints - CritTraj{j}(end)));
    if delta_s<1e-2 && delta_f<1e-2 
        Adjacency(start,finish)=j;
    else
        Adjacency(start,end)=j;
    end
end
CritTraj=CritTraj(1:countcrit); 
end