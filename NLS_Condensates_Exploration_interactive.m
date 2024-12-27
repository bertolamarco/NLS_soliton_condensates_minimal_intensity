
% The data that we want shared between callbacks is in setappdata(gca), it
% contains all the arguments of the main computatio function (Chebotarov or
% Critical energies).
% The external field is encoded in phi, representing a polynomial without
% constant term and real coefficients. By default phi=[1,0], representing
% the external field z. The script creates 4 random anchor points in the
% upper half plane. They can be moved around with the mouse. The buttons at
% the bottom of the figure allow to change the number of stagnation points
% (L=L+1, L=L-1) respectively, and to add or remove anchor points (N=N+1,
% N=N-1), respectively. 



% Copyright Marco Bertola 2024.
% marco.bertola@concordia.ca

N=4;
E = randn(1,N) + 1i*abs(randn(1,N));
fig = sfigure(1);
clf;
hold on;

L=1;
sroot=[];
droot=[];
epsilon=0.001;
Length=3;
phi=[1,0]; %External field is z;


bplay = uicontrol(fig, ...
    "String","Recompute", ...
    "Position",[310 10 80 20]...
    ,'Callback',@(src,event) Recompute(src,event, fig), 'Interruptible','off');

bplus = uicontrol(fig, ...
    "String","L=L+1", ...
    "Position",[70 10 50 20]...
    ,'Callback',@(src,event) ChangeL(src,event,1,fig));

bminus = uicontrol(fig, ...
    "String","L=L-1", ...
     "Position",[130 10 50 20]...
    ,'Callback',@(src,event) ChangeL(src,event,-1,fig));

bplus = uicontrol(fig, ...
    "String","N=N+1", ...
    "Position",[190 10 50 20]...
    ,'Callback',@(src,event) ChangeN(src,event,1,fig));

bminus = uicontrol(fig, ...
    "String","N=N-1", ...
     "Position",[250 10 50 20]...
    ,'Callback',@(src,event) ChangeN(src,event,-1,fig));




AnchorHandles = zeros(size(E));
for j=1: size(E,2)
    AnchorHandles(j) = plot(E(j),'.r','MarkerSize',22, 'ButtonDownFcn',  @(src,event) Grab(src,event));
    setappdata(AnchorHandles(j),'Identity',j);
    hold on;
end
[sroot,droot,Amp,Bcuts, Objects] = Critical_Energies_NLS(phi,L,E,true,sroot,droot,epsilon, Length);
L = size(sroot,2);
setappdata(gca,'AnchorSet',E);
setappdata(gca,'AnchorHandles',AnchorHandles);
setappdata(gca,'NumStagPts',L);
setappdata(gca,'NumAnchors',N);
setappdata(gca,'InitDroot',droot);
setappdata(gca,'InitSroot',sroot);
setappdata(gca,'epsilon',epsilon);
setappdata(gca,'Length',Length);
setappdata(gca,'phi',phi);
%Objects = Chebotarov(E,L);
title({['The Dirichlet Energy is ',num2str(Amp)]; ['We have L=',num2str(L),' stagnation points'] },'FontSize',18);
setappdata(gca,'Objects',Objects);
%plot([0,0.4i,0.4+0.4i], 'ButtonDownFcn',  @(src,event) Grab(src,event));
grid on;
set(gca, 'XLimMode','manual','YLimMode','manual');
set(gca,'Xlim',[-4,4]);
set(gca,'Ylim',[0,4]);
arrayfun(@(p) uistack(p,'top'), getappdata(gca, 'AnchorHandles')); %Set the anchors to the top of the graphics stack so that we can still mouse click on them!


%%% END OF SCRIPT %%%



function Recompute(src,event, fig)
    if ~isvalid(fig)
        delete(src.Parent);
    else
        ax =  findobj(fig.Children, 'Type', 'axes'); % Should be the only axes handle in the figure!
        E = getappdata(ax,'AnchorSet');
        L = getappdata(ax,'NumStagPts');
        disp('Do Something but Nothing');
        oldObj=getappdata(ax, 'Objects');
        delete(oldObj);
        sroot=getappdata(ax,'InitSroot');
        droot=getappdata(ax,'InitDroot');
        epsilon=getappdata(ax,'epsilon');
        Length=getappdata(ax,'Length');
        phi=getappdata(ax,'phi');
        
        %Objects = Chebotarov(E,L);
        [sroot,droot,Amp,Bcuts, Objects] = Critical_Energies_NLS(phi,L,E,true,sroot,droot,epsilon, Length);
        L=size(sroot,2);
        title({['The Dirichlet Energy is ',num2str(Amp)]; ['We have L=',num2str(L),' stagnation points'] },'FontSize',22);

        setappdata(gca,'NumStagPts',L);
        setappdata(gca,'InitDroot',droot);
        setappdata(gca,'InitSroot',sroot);
        setappdata(gca,'epsilon',epsilon);
        setappdata(gca,'Length',Length);
        setappdata(ax,'Objects',Objects);
        arrayfun(@(p) uistack(p,'top'), getappdata(ax, 'AnchorHandles')); %Set the anchors to the top of the graphics stack so that we can still mouse click on them!
    end
end

function ChangeL(src,event, val, fig)
    if ~isvalid(fig)
        delete(src.Parent);
    else
        ax =  findobj(fig.Children, 'Type', 'axes'); % Should be the only axes handle in the figure!
        L= getappdata(ax,'NumStagPts');
        L=max(L+val,0);
        setappdata(ax,'NumStagPts',L);
        disp(['L=',num2str(L)]);
        %Since we changed the number of stagnation points, the next
        %"Recompute" will have to run from scratch, so we empty sroot,
        %droot.
        setappdata(gca,'InitDroot',[]);
        setappdata(gca,'InitSroot',[]);
        title (['L changed to ',num2str(L),'. Click Recompute']);
    end

end

function ChangeN(src,event, val, fig)
    if ~isvalid(fig)
        delete(src.Parent);
    else
        ax =  findobj(fig.Children, 'Type', 'axes'); % Should be the only axes handle in the figure!
        N= getappdata(ax,'NumAnchors');
        N=N+val;
        E = getappdata(ax,'AnchorSet');
        AnchorHandles=getappdata(gca,'AnchorHandles');

        if val==-1
            E = E(1:end-1);
            delete(AnchorHandles(end));
            AnchorHandles=AnchorHandles(1:end-1);
        else
            E =[E, randn(1) + 1i*abs(randn(1))];
            AnchorHandles=[AnchorHandles,0];
            AnchorHandles(end) = plot(E(end),'.r','MarkerSize',22, 'ButtonDownFcn',  @(src,event) Grab(src,event));
            setappdata(AnchorHandles(end),'Identity',size(E,2));
        end

        setappdata(ax,'NumAnchors',N);
        setappdata(ax,'AnchorHandles',AnchorHandles);
        setappdata(ax,'AnchorSet',E);
        disp(['N=',num2str(N)]);
        %Since we changed the number of stagnation points, the next
        %"Recompute" will have to run from scratch, so we empty sroot,
        %droot.
        setappdata(gca,'InitDroot',[]);
        setappdata(gca,'InitSroot',[]);
        title (['N changed to ',num2str(N),'. Move the anchors and click Recompute']);
    end

end



function Grab(src,event)
% If we click on an object we set the "MotionFnc" to move it around with Drag, until
% we release the click and reset the MotionFnc with Release
parfig=gcbf;
    if src.Type == 'line'
        obj=src;
        setappdata(obj,'OldButtonUpFcn', get(parfig,'WindowButtonUpFcn'));
        setappdata(obj,'GrabPoint',gca().CurrentPoint);
        %setappdata(obj,'OldButtonDownFcn', get(parfig,'WindowButtonDownFcn'));
        setappdata(obj,'OldMotionFcn', get(parfig,'WindowButtonMotionFcn'));
        set(parfig, 'WindowButtonMotionFcn', @(src,event) Drag(src,event,obj) );
        set(parfig, 'WindowButtonUpFcn', @(src,event) Release(src,event,obj));
    end

end


function Drag(src,event,obj)
    ax=gca;
    OldPt = getappdata(obj,'GrabPoint');
    NewPt = ax.CurrentPoint;
    set(obj,'XData', get(obj,'XData') + NewPt(1,1) -  OldPt(1,1));
    set(obj,'YData', get(obj,'YData') + NewPt(1,2) -  OldPt(1,2));
    setappdata(obj,'GrabPoint', NewPt);
    Id = getappdata(obj,'Identity');
    E = getappdata(ax,'AnchorSet');
    E(Id) = NewPt(1,1)+1i*NewPt(1,2);
    setappdata(ax,'AnchorSet',E);
end

function Release(src,event,obj)
    set(src,'WindowButtonUpFcn', getappdata(src,'OldButtonUpFcn'));
    %set(obj,'ButtonDownFcn',getappdata(src,'OldButtonDownFcn'));
    set(src,'WindowButtonMotionFcn',getappdata(src,'OldMotionFcn'));

    ax=gca;
    % Make sure that the AnchorSet has the same coords as the one displayed
    AnchorHandles = getappdata(ax,'AnchorHandles');
    AnchorSet = getappdata(ax,'AnchorSet');

    for j=1:size(AnchorHandles,2)
        set(AnchorHandles(j),'XData', real(AnchorSet(j)));
        set(AnchorHandles(j),'YData', imag(AnchorSet(j)));
    end
     
    
    
end

function displayCoordinates(src,~,ax)
src.MarkerEdgeColor = rand(1,3);
disp(ax.CurrentPoint(1,1:2))
end