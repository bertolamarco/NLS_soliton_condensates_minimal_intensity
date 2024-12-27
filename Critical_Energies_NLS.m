
function [sroot,droot,Amp, BCuts, Objects]=Critical_Energies_NLS(phi,L,e, Draw,sroot,droot,epsilon,Length)
%% [sroot,droot,Amp,BCuts, Objects]=Critical_Measures_NLS(phi, L,e, Draw,sroot, droot, epsilon, Length)
%% Finds the double and simple zeroes of the generalized quasienergy quadratic differential (dp)^2 
%% corresponding to the critical measure with endpoints e=(e_1,...e_N) and external field given by
%% the real-coefficient polynomial phi (without constant coefficient!) of degree $R$
%% in the upper half plane. 
%% The second optional argument L is the number of simple zeroes (L<=N)
%% The output is sroot= double roots, droot=simple roots, Amp=Amplitude of the corresponding critical measure, namely
%% the dp = phi'(z) -Amp/z^2 + ...
%% The form of dp is dp = S(z) \sqrt(D(z))/\sqrt(A(z)) 
%% where $D(z) = \prod_{j=1}^{N+R-1-L} (z-d_j)(z-\overline d_j)$, 
%% $S(z) =\prod_{j=1}^{L} (z-s_j)$
%% $A(z) = \prod_{j=1}^{N} (z-e_j)(z-\overline e_j)$
%% Draw is a flag to decide if we want to draw the figures
%% sroot are the roots of S, droot are the roots of D in the upper half plane.
%% Amp is the amplitude of the NLS gas.
%% BCuts should be the branch-cuts going along the critical trajectories.
%% Objects is a list of handles to all graphical objects plotted, in case we want to manipulate them externally.
%% Copyright Marco Bertola 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Initialization section: %%%%%%%%%%
    Precision =1e-8;
    MAXITER=300; %Seems enough.
    if ~exist ('phi','var')
        %phi = [1 , 0 ,0]; % Phi = z^2
        phi = [1,-2,0];
    end
    phi(end)=0; %we make sure that the constant coefficients is zero!

    if ~exist ('epsilon','var')
        epsilon =0.003;
    end
    if ~exist ('Length','var')
        Length =4;
    end
  
    if ~exist ('e','var') 
        %e =0.37*[-2,-1,1,2,3,7] + (2+[1,1.9,0.94,1.2,1,-1])*1i; %This violates the conjecture!
        %e=[-1,0.8,1]*2 + [3,4,3]*1i/4;
        e = [-2, -1,0,1,2]/2-1/4 + 1i*[0.3, 1,2,2,1];
        %e = [-2 2] +2i*[0.4,0.5];
        %e=[-1,1,0.6] + [2,2,1]*1i;
        %e=[ 1i+0.5 ,1i - 1.3];
    end
    e= sort(e,'ComparisonMethod','real');
    N= size(e,2);
    R = size(phi,2)-1; % the degree of the external potential

    if ~exist('L','var')
        L=1;
        %By default we look simply for the standard ZS critical measure
        %with the maximum number of stagnation points. In general we may
        %have L stagnation points and d = N+R-L additional branchpoints in
        %the upper half plane. 
    end
    if ~exist('Draw','var')
        Draw=1; % The default is to draw.
    end
    if L<0 
        L=0;
    end
    if size(e(imag(e)<=0),2)>0
        warning('The endpoints must be in the upper half plane!');
        e = real(e) + 1i*abs(imag(e));
    end
    if L>N+R-1
        warning('The value of L should be no greater than the number of endpoints');
        L=N;
    end
    numzero=1e-3;
    ITER = 0;
    Adaptive=1;
    %% Change Adaptive to a bigger value if it does not converge!
    e=sort(e,'ComparisonMethod','real');
    if exist('droot','var')
        if isempty(droot)
            clear('droot');
        end
    end
    if exist('sroot','var')
        if isempty(sroot)
            clear('sroot');
        end
    end

    M = N+R-1-L;%This is the number of simple roots in D(z)
    if ~exist('droot','var')
        % We choose the roots of D randomly in the convex hull of E and the
        % real axis. 
       
        pp=cat(2,[min(real(e)), max(real(e))], e);
        s0 = zeros(M,N+2);
        droot=zeros(1,M);

        for j=1:M
            MM = rand(1,N+2);
            s0(j,:) = pp.*MM/sum(MM);
            droot(j) = sum(pp.*MM/sum(MM));
        end
    end
    if ~exist('sroot','var') % We choose the roots of S randomly along the real axis at first.
        s0=2*rand(1,N+R-1-M);
        sroot =min(real(e))*s0 + (1-s0)*max(real(e)) ;
    end
    if size(sroot,2)~=L || size(droot,2)~= M
        disp('The number of sroot or droot is incompatible with the passed L. Fixing myself...');
        L= size(sroot,2);
        M=N+R-1 -L;
    end
    droot= sort(droot,'ComparisonMethod','real');

   
    %%%%%%%%%%%%%%%%%%%%% End of Initialization section %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    D=poly(cat(2,droot ,conj(droot)));
    Fun=2;
    Funold=Fun;
    Bounce=0; %This is to detect flip of sheet and wrong trend;
    Sign=1;
    T=polyder(phi);
    S = T(1)^2*poly(sroot);

    if Draw
        dd=plot(droot,'.k','MarkerSize',10);
        hold on;
        ee=plot(e,'.r','MarkerSize',20);
        ss=plot(real(sroot),imag(sroot), '.g', 'MarkerSize',12);
    end
    axis ([min(real(e))-1,max(real(e))+1,0,max(imag(e))+1]); hold on;
    %end initialization;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%MAIN LOOP %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while(ITER<30 || (Fun>Precision && ITER<MAXITER)) 
        %Main loop
        [Inty0, P]  = Cycles(T(1), e,sroot,droot);
        Inty=Inty0;
        %g=N+size(droot,2)-1;
        tmp1= cat(2,T(2:end),0);
        Inty(1:R) = (Inty(1:R)- 2*pi*tmp1') ;
        %\dot t_\ell = -(t_ell-phi_ell); ell=1,..,R-1 plus t_0 which is the
        %residue.
        %The leading coefficients should be zero (numerically) because we
        %already fixed the leading coefficient.

        Fun =sqrt( sum(Inty.^2));
        if Funold<Fun
            Bounce=Bounce+1;
        else
            Bounce=0;
        end
        if Bounce>3
            Sign=-Sign;
            Bounce=0;
        end
        Funold=Fun;
        %we have y = S sqrt(D)/sqrt(A) so that
        %doty = (dot S D +1/2*S dotD)/sqrt(AD)
        %This is a holomorphic differential dotF/sqrt(AD) such that the
        %imaginary periods are exactly minus Inty.
        % We call dotF the numerator of this expression: it is a polynomial
        % of degree 2M+L-1
        dotF = -Sign*real(P^(-1)*Inty);
        %dotF=dotF(1:end-1).';
        %%%  DEGREE of dotF should be g+R-1
        %dotF= dotF(end:-1:1);
        %dotF=cat(1,0,dotF).';
        dotF = dotF .'; %should be a row vector for the computation below using the Sylvester matrix.
        % We have to solve the polynomial equation dotS*D + S*dotD/2 = dotF for
        % dotS, dotD, given S,D,dotF. note that the degree of dotS<=degree S-1
        % (similarly for dotD and D). This equation leads to a linear system
        % for the coefficients of dotS, dotD whose matrix is the Sylvester
        % matrix.
        MM=mysylvester(S,D);
        V=dotF*MM^(-1);
        dotS = V(2*M+1:end);
        dotD = 2*V(1:2*M);
        %padding to correct size
        dotD= cat(2,0,dotD);
        dotS= cat(2,0,dotS);   
        % if sum(abs(conv(S,dotD)/2 + conv(D,dotS)- dotF))>1e-9
              % error('Cretino!');
        % end
        c1=sum(abs(dotS))/sum(abs(S));
        c2=sum(abs(dotD))/sum(abs(D));
        dt=(1/2+( c1/(1+c1) + c2/(1+c2))/4)/Adaptive; %Choice of step. 
    
        if(size(D,2)>1)
            D=D+dt*dotD;
        end
        S=S+dt*dotS;
    
        %We need to check if there are new real roots from D
        droot=roots(D); %row vector, please...
        droot= sort(droot,'ComparisonMethod','real');
        droot=droot.';
    
        d_up = droot(imag(droot)>numzero);
        d_R=real(droot(imag(droot)<=numzero & imag(droot)>=0));
    
        %We need to check that we did not have a bounce of complex roots on the
        %real axis; in this case we have to force-merge
        if ((size(D,2)-1)/2 - size(d_up,2))< size(d_R,2)
            d_R = sort(d_R);
            d_R = (d_R(1:2:end)+d_R(2:2:end))/2;
        end
        %d_R = (d_R(1:2:end)+d_R(2:2:end))/2;
        S = conv(S,poly(d_R));
        %Lold=L;
        %L=size(S,2)-1;
        sroot = roots(S);
        sroot = sort(sroot,'ComparisonMethod','real');
        sroot = sroot.';
    
        %real zeros of D must be double!
        droot=d_up; %Only the roots that did not merge with the conjugate on the real axis remain in D; the other go to S.
        M= size(droot,2);
        ITER = ITER+1;
        if ITER>50 && mod(ITER,10)==0
            Adaptive=min(Adaptive+1,1);
        end

        % We should check also if two zeros of D have merged in the upper
        % half plane because this causes instability;

        dist = abs( ones(M,1)*droot - droot.'*ones (1,M)) + eye(M);
        
        if min(min(dist))<numzero
            [J,K]= find(dist<numzero); %J is of even length because dist is symmetric!
            J = J(1:end/2);
            %K = K(1:end/2);

            S = conv(S,poly(cat(2,droot(J),conj(droot(J)) )));
            droot(K)=[];
            M = size(droot,2);
            %Lold=L;
            %L=size(S,2)-1;
            sroot = roots(S);
            sroot = sort(sroot,'ComparisonMethod','real');
            sroot = sroot.';

        end
        D = poly(cat(2,droot,conj(droot)));
    
        %------------------------
        %Drawing section
        if Draw 
            sroot=roots(S).';
            delete(dd);
            delete(ee);
            delete(ss);
            ss=plot(real(sroot),imag(sroot), '.g','MarkerSize',10);hold on;
            dd=plot(real(cat(2,droot,conj(droot))),imag(cat(2,droot,conj(droot))) ,'.k','MarkerSize',10);hold on;
            ee=plot(cat(2,e,conj(e)),'.r','MarkerSize',10);

            title(['Fun=',num2str(Fun), '. Iteration=', num2str(ITER)]);
            drawnow;
        else
            disp(['L^2 norm of imaginary periods=',num2str(Fun),'. Iteration=',num2str(ITER)]);
        end
        %------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% End Main Loop %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Clean up and exit %%%%
    %disp(['Converged in ',num2str(ITER),' iterations'])
    %%%% We now plot the critical trajectories.
    if(ITER==MAXITER)
        disp(['Max iter! Tolerance=',num2str(Fun)]);
        if Fun>1
            stop;
        end
    end
    
    y=Qmom(T(1),e,sroot,droot);
    W=max(abs(e))*2;
    Amp=real(-integral( @(z) polyval(phi,z).*y(z), W,W,'Waypoints',[W*1i,-W,-W*1i])/(2i*pi));

    %Mass= real(integral( @(z) y(z), -W,-W,'Waypoints',[W,W+1i*W,-W+W*1i])); %It should be numerically real
    %disp(['The mass is =',num2str(real(Mass))]);        
    disp(['The energy is =',num2str(Amp)]);
    %Now the amplitude of the ZS spectrum.
    %[~,AmpZS,~]=QuasimomentumNLS(e,0);
    %The amplitude needs to be computed as an integral, not as a limit, because
    %of numerical stability.    
    
    %%% Plot all the scenarios and done!
    if Draw
        SSS = cat(2,sroot,droot, conj(droot));
        EEE = cat(2,e,conj(e));
        [CTraj,Adjacency,TypeCritpts, ~]=CriticalTrajectories_RationalNODRAW (1,SSS ,cat(2,2*ones(size(sroot)), ones(size(cat(2,droot,conj(droot))))) ,...
            EEE,ones(size(EEE)),epsilon,Length);
        Traj=zeros(size(CTraj));
        color=zeros(size(CTraj));

        for j=1: size(CTraj,2)
            [start,~]=find(Adjacency==j); % The j trajectory connects the critical point l to k.
            color(j)='g';
            if ~isempty(start) 
                if rem(TypeCritpts(start),2)==1 
                    color(j) = 'b';             
                else
                    color(j)= 'g';
                end
            end
        end
        Traj(color=='g') = cellfun(@(p) plot(real(p),imag(p),'g','LineWidth',1),CTraj(color=='g'));
        Traj(color=='b') = cellfun(@(p) plot(real(p),imag(p),'b','LineWidth',3),CTraj(color=='b'));

        % now we find the trajectories connecting e's or an e to the real
        % axis;

        BCuts=CTraj;
        count=0;
        Taken= zeros(2,size(BCuts,2));

        for j=1:size(CTraj,2)
            gamma=CTraj{j};
            [d,jj]=min(abs(gamma(1)-e));
            if d<1e-2 % it starts  at e(jj)
                %Check if it ends at another e
                [d_d,kk]=min(abs(gamma(end)-e));
                if d_d<1e-2  % it does; but was it already considered?
                    if isempty(Taken(Taken(1,1:count)==kk))...
                            && isempty(Taken(Taken(2,1:count)==kk))%it ends at e(kk) and it wasn't counted before
                        count=count+1;
                        BCuts{count} = gamma;
                        Taken(:,count)=[jj;kk]; %from e(jj) to e(kk);
                    end
                else
                    %It may be a spurious count that goes to a saddle
                    %point; let's check
                    if imag(gamma(end))<1e-2
                        %It goes to the real axis. We find which sroot it goes to and truncate and store.
                        if isempty(Taken(Taken(1,1:count)==jj))...
                                && isempty(Taken(Taken(2,1:count)==jj))
                            count=count+1;
                            [ddd,aaa]=min ((abs((ones(size(sroot,2),1)*gamma - sroot.'*ones(1,size(gamma,2)) ))),[],2);
                            [~,aa]=min(ddd);
                            BCuts{count} = [gamma(1:aaa(aa)), sroot(aa)];
                            Taken(:,count)=[jj;0]; % From e(jj) to the real axis.
                        end
                    end
                end
            end

        end
        BCuts=BCuts(1:count);

        hold on;grid on;
        delete(ee);
        delete(dd);
        delete(ss);
        ee=plot(cat(2,e,conj(e)),'.r','MarkerSize',10);
        dd=plot(cat(2,droot,conj(droot)),'.k','MarkerSize',12);
        SSS=cat(2,sroot,conj(sroot));
        ss= plot(real(SSS), imag(SSS), '.g','MarkerSize',12);
        win = [ min(real(e))-2, max(real(e))+2,...
            0,...%-max(imag(e))-3,
            max(imag(e))+2];
        axis (win);
        drawnow;
        Objects=[Traj,dd,ss,ee];

       
          axis equal
          axis(win);
          name = 'phi=';
          for j =1:size(phi,2)-1
              name= cat(2, name ,[num2str(phi(j)),' z^',num2str(R+1-j),'+']);
          end
          name = name(1:end-1);
          name = cat(2,name,'.pdf');

          title('');
          %exportgraphics(gcf,name,'ContentType','vector');
          disp(['droot=',num2str(droot)]);
          disp(['sroot=',num2str(sroot)]);
          
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% END OF SCRIPT %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function y=Qmom(lcoeff,e,sroot,droot)
    numzero=1e-3;
    N=size(e,2);
    M=size(droot(imag(droot)>numzero),2);
    d_up=droot(imag(droot)>numzero);
    if size(d_up)~= size(droot)
        error('Ough');
    end
    t1= cat(2,e,d_up);
    t2= cat(2,conj(e), conj(d_up));
    %Vertical cuts;
    
    Cuts= cat(2,t1.', t2.');
    
    t1 = cat(2,d_up,conj(d_up),sroot);
    t1=real(poly(t1)); %it should be real, but numeric fluctuations may give small imaginary part, which we kill...
    T = @(z) lcoeff*polyval(t1,z);

    % R is the radical of the Riemann surface; y = T/R dz
    R = @(z) (1i)^(N+M)*...
        prod(...
        sqrt( -1i*(ones(N+M,1)*z-Cuts(:,1)*ones(size(z))) )...
        .*sqrt( -1i*(ones(N+M,1)*z-Cuts(:,2)*ones(size(z))) ));
    if (N+M)==1 %it means there is only one e and no droot
        R= @(z) 1i*sqrt(-1i*(z-e(1))).*sqrt(-1i*(z-conj(e)));
    end

    y= @(z) T(z)./R(z);
end

function [Inty,PowerPeriods] = Cycles (tR, e,sroot,droot)
%Apparently, due to Schwartz symmetry,
% Computes the integral of y dz  along the two cycles as well as the holomorphic ones
%%It is assumed that the odd d,e's are in the upper half plane, ordered by
%%real part.
%sroot are the double roots, droot are the simple roots (not real)
numzero=1e-3;
%tR is the leading coefficient of Phi

N=size(e,2);
M=size(droot(imag(droot)>numzero),2);
d_up=droot(imag(droot)>numzero); %the roots of D in the upper half plane
if size(d_up)~= size(droot)
    error('Ough');
end
t1= cat(2,e,d_up);
t2= cat(2,conj(e), conj(d_up)); % the roots of E(z)D(z) (2N+2M of them). This is the Riemann surface, g= N+M-1
%Vertical cuts;

Cuts= cat(2,t1.', t2.');

t1 = cat(2,d_up,conj(d_up),sroot);
t1=poly(t1);
T = @(z) tR*polyval(t1,z); % the polynomial S D

% Rad is the radical of the Riemann surface; y = T/R dz

Rad = @(z) (1i)^(N+M)* prod(...
    sqrt( -1i*(ones(N+M,1)*z-Cuts(:,1)*ones(size(z))) )...
    .*sqrt( -1i*(ones(N+M,1)*z-Cuts(:,2)*ones(size(z)))) ) ;
%end
if (N+M)==1 %it means there is only one e and no droot
    Rad= @(z) 1i*sqrt(-1i*(z-e)).*sqrt(-1i*(z-conj(e)));
end
y= @(z)  T(z)./Rad(z);

r=0.0007;
g= N+M-1;
R = size(sroot,2)+ M+1-N; % the order of the pole (-1) and degree of phi.

Path=cell(g+R,1);
% each path has a  component on the first and on the second sheet (possibly none of the latter)
%Don't put any endpoint on a cut!
for j=1:g
    Path{j+R}   =[Cuts(j,1)+r*exp(1i*pi*( (0:1/5:1))), Cuts(j,2)+r*exp(1i*pi*( (1:1/5:2))), Cuts(j,1)+r] ;
    %A-cycles
end
KK=30; %something large for the circle around infinity
for j=1:R
    Path{j}   =KK*[1+1i,-1+1i,-1-1i,1-1i, 1+1i] ; % circle around infinity, repeated R+1 times
end

PowerPeriods= zeros(g+R,g+R); 
%This will store the (purely imaginary, thanks to Schwartz symmetry)
%A-periods around the vertical bands. 
Inty=zeros(g+R,1); %


for k=1:g+R  
    F = @(z) z.^(g+R-k)./Rad(z);
    for j=1:R %% added -1 to g+R
        PowerPeriods(j, k) =  Hyperelliptic_Integral (@(z) F(z)./z.^(R-j),Path{1}, Cuts); 
        %The first R+1 paths are all circles. The R+1 first minor is made
        %of the residues of second-kind differentials at infty.
    end
    for j=R+1:g+R
        PowerPeriods(j, k) =  Hyperelliptic_Integral (@(z) F(z),Path{j}, Cuts);
    end
end

for j=1:R
    Inty(j) =  Hyperelliptic_Integral (@(z) y(z)./z.^(R-j),Path{1}, Cuts);
end
for j=R+1:g+R
    Inty(j) =  Hyperelliptic_Integral (@(z) y(z),Path{j}, Cuts);
end
if( sum(sum(real(PowerPeriods).^2))>1e-5)
    error('There are some real part too big here....!');
end
PowerPeriods=imag(PowerPeriods);
Inty = imag(Inty);
end

function Syl= mysylvester(S,D)
 % Given the two polynomials S, D it returns the Sylvester matrix. No
 % sanity check is done. YOLO! NOTE: if S,D are constants, it returns a 1x1
 % matrix and not an empty matrix. But why would you compute the Sylvester
 % matrix of two constants?

 s = size(S,2);
 d = size(D,2);
 
 %Syl = zeros(s+d-2,s+d-2);%
 if s==1
     Syl=toeplitz([S(1);zeros(d-2,1)], [S,zeros(1,d-2)]);
 else
     if d==1
         Syl=toeplitz([D(1);zeros(s-2,1)], [D,zeros(1,s-2)]);
     else
         Syl=[toeplitz([S(1);zeros(d-2,1)], [S,zeros(1,d-2)]);
             toeplitz([D(1);zeros(s-2,1)], [D,zeros(1,s-2)])];
     end
 end

end