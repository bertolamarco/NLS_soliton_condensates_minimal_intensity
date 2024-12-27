function J=Hyperelliptic_Integral(F, Path,Cuts,varargin)
%
% August 2024
% (c) Marco Bertola
% email: Marco.Bertola@concordia.ca
%% J=Hyperelliptic_Integral(F, Path,Cuts) computes the integral of F along Path,
%% where F is a function on the hyperelliptic Riemann surface obtained by gluing 
%% two copies of C along Cuts. It is assumed that F is odd  on this hyperelliptic surface and no check is done.
%
silent=0;
if ~isempty(varargin)
    if varargin{1}=="Silent" && varargin{2}==true
        silent=1;
    end
end
PathSheet = SheetPartition(Path,Cuts);

if PathSheet{1}{1}(1)~=PathSheet{1}{end}(2) && ~silent
    disp('Warning: the path does not close on the Riemann surface...');
end
J=0;
for pp= PathSheet{1} %pieces on the first sheet
    J = J+ integral( @(z)  F(z), pp{1}(1), pp{1}(end), 'Waypoints',pp{1});
end
for pp=PathSheet{2} %pieces on the second sheet
    J = J- integral( @(z)  F(z), pp{1}(1), pp{1}(end), 'Waypoints',pp{1});
end
end

function Paths=SheetPartition(Q, Cuts)
% Splits the polygonal path Q into subsegments according to the sheets; each time we cross a cut we jump on other sheet.
%Cuts must be a Nx2 array (starting and ending point)
%Q must be a vector specifying the nodes of the path.
%It returns a 2x1 cell where Paths{sheet} (sheet=1,2) is in turn a cell
%containing all segments on sheet "sheet". 
Ncuts = size(Cuts,1);
Paths=cell(2,1); %These will contain the segments on the first and second sheet, respectively
sheet=1; %We start on sheet 1
 %The initial point is on the first sheet; 
for h=1:size(Q,2)-1
   Flags=zeros(Ncuts,3); 
   for j=1:Ncuts
       [flag,pt,s]= Intersect(Q(h), Q(h+1), Cuts(j,1),Cuts(j,2));
       if flag==1
           Flags(j,:)=[flag,pt,s];
       end %Now Flags contains all intersections of the segment [Q(h),Q(h+1)] with all cuts; 
   end
   [~,J]=sort(Flags(:,3));
   Flags= Flags(J,:);
   %Prune it down to only intersections
   Intpts= Flags(Flags(:,1)==1,2);
   L=[Q(h),Intpts.',Q(h+1)];
   %The segment crosses at least one cut.
  
   for g=1:size(Intpts,1)
       Paths{sheet}{end+1}=[L(g), L(g+1)];
       sheet=3-sheet; %Flip sheet;
   end
   if(size(Intpts,1)==0) 
       Paths{sheet}{end+1}=[Q(h), Q(h+1)];
   else
       Paths{sheet}{end+1}=[L(g+1), L(g+2)];
   end
   %Hopefully the last or initial  point is not exactly on a cut!
end
end

function [flag,pt,s1] = Intersect(z1,z2,w1,w2)
%Checks if the segments [z1,z2] and [w1,w2] intersect; they do iff s1 and
%s2 are both in [0,1].
flag=0;pt=0;
Delta = imag( (z1-z2)*conj(w1-w2));
s1= 1-imag( (w2-z2)*conj(w1-w2))/Delta;
s2= -imag( conj(w2-z2)*(z1-z2))/Delta;
if ( s1>=0 && s1<=1 && s2>=0 &&s2<=1)
    flag=1;
    pt= s1*z2+(1-s1)*z1;
end
end
