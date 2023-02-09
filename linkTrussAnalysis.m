close all
clear
clc
%% Specify Tiles (corner locations)
% Allows any polygon shapes
% Planar tiling
% tiles(1) = struct('nodes',[0 0 0; 1 0 0; 1 1 0; 0 1 0],'nodeInds',[],'ndir',[]);
% tiles(2) = struct('nodes',[1 1 0; 2 1 0; 2 2 0; 1 2 0],'nodeInds',[],'ndir',[]);
% tiles(3) = struct('nodes',[0 2 0; 1 2 0; 1 3 0; 0 3 0],'nodeInds',[],'ndir',[]);
% tiles(4) = struct('nodes',[-1 1 0; 0 1 0; 0 2 0; -1 2 0],'nodeInds',[],'ndir',[]);
% tiles(4) = struct('nodes',[-1 1 0; 0 1 0; 0 2 0; -1 2 0; -1.5 1.25 0],'nodeInds',[],'ndir',[]); % pentagon
% Tessellation: unit cell: 2x2.
% marr = 3;
% narr = 3;
% a = 1;
% b = 2;
% tiles(2*marr*narr) = struct('nodes',[],'nodeInds',[],'ndir',[]);
% baseNodes = [0 0 0; a 0 0; a b 0; 0 b 0];
% tiles(1) = struct('nodes',baseNodes,'nodeInds',[],'ndir',[]);
% tiles(2) = struct('nodes',baseNodes+[1 0 0]*a + [0 1 0]*b,'nodeInds',[],'ndir',[]);
% tileCount = 3;
% xRepEachSide = floor((marr-1)/2);
% yRepEachSide = floor((narr-1)/2);
% for i = -xRepEachSide:(xRepEachSide+1-mod(marr,2))
%     for j = -yRepEachSide:(yRepEachSide+1-mod(narr,2))
%         if ~( i == 0 && j == 0)
%             xshift = 2*a*i;
%             yshift = 2*b*j;
%             tiles(tileCount) = struct('nodes',baseNodes+[1 0 0]*xshift + [0 1 0]*yshift,'nodeInds',[],'ndir',[]);
%             tiles(tileCount+1) = struct('nodes',baseNodes+[1 0 0]*(xshift+a) + [0 1 0]*(yshift+b),'nodeInds',[],'ndir',[]);
%             tileCount = tileCount + 2;
%         end
%     end
% end
% Tessellation with 1DOF: unit cell: 4x4 modified.
% marr = 2;
% narr = 2;
% a = 1;
% b = 1;
% tiles(8*marr*narr) = struct('nodes',[],'nodeInds',[],'ndir',[]);
% baseNodes = [0 0 0; a 0 0; a b 0; 0 b 0];
% tiles(1) = struct('nodes',baseNodes,'nodeInds',[],'ndir',[]);
% tiles(2) = struct('nodes',baseNodes-[1 0 0]*a+[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(3) = struct('nodes',baseNodes+[1 0 0]*a+[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(4) = struct('nodes',baseNodes+2*[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(5) = struct('nodes',baseNodes+2*[1 0 0]*a+2*[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(6) = struct('nodes',baseNodes+[1 0 0]*a+3*[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(7) = struct('nodes',baseNodes+3*[1 0 0]*a+3*[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tiles(8) = struct('nodes',baseNodes+2*[1 0 0]*a+4*[0 1 0]*b,'nodeInds',[],'ndir',[]);
% tileCount = 9;
% xRepEachSide = floor((marr-1)/2);
% yRepEachSide = floor((narr-1)/2);
% for i = -xRepEachSide:(xRepEachSide+1-mod(marr,2))
%     for j = -yRepEachSide:(yRepEachSide+1-mod(narr,2))
%         if ~( i == 0 && j == 0)
%             xshift = 4*a*i;
%             yshift = 4*b*j;
% tiles(tileCount) = struct('nodes',baseNodes+[1 0 0]*(xshift)+[0 1 0]*(yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+1) = struct('nodes',baseNodes-[1 0 0]*(a-xshift)+[0 1 0]*(b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+2) = struct('nodes',baseNodes+[1 0 0]*(a+xshift)+[0 1 0]*(b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+3) = struct('nodes',baseNodes+[1 0 0]*xshift+[0 1 0]*(2*b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+4) = struct('nodes',baseNodes+[1 0 0]*(2*a+xshift)+[0 1 0]*(2*b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+5) = struct('nodes',baseNodes+[1 0 0]*(a+xshift)+[0 1 0]*(3*b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+6) = struct('nodes',baseNodes+[1 0 0]*(3*a+xshift)+[0 1 0]*(3*b+yshift),'nodeInds',[],'ndir',[]);
% tiles(tileCount+7) = struct('nodes',baseNodes+[1 0 0]*(2*a+xshift)+[0 1 0]*(4*b+yshift),'nodeInds',[],'ndir',[]);
% %             tiles(tileCount) = struct('nodes',baseNodes+[1 0 0]*xshift + [0 1 0]*yshift,'nodeInds',[],'ndir',[]);
% %             tiles(tileCount+1) = struct('nodes',baseNodes+[1 0 0]*(xshift+a) + [0 1 0]*(yshift+b),'nodeInds',[],'ndir',[]);
%             tileCount = tileCount + 8;
%         end
%     end
% end

% Singly curved tiling
theta = deg2rad(30);
tile_pre_Rotate = [0 0 0;
    1 -1 0;
    2 0 0;
    1 1 0;
    2 0 0;
    3 -1 0;
    4 0 0;
    3 1 0;
    0 1+cos(theta) -sin(theta);
    1 1 0;
    2 1+cos(theta) -sin(theta);
    1 1+2*cos(theta) -2*sin(theta);
    2 1+cos(theta) -sin(theta);
    3 1 0;
    4 1+cos(theta) -sin(theta);
    3 1+2*cos(theta) -2*sin(theta)];
rotAngDeg = 45;
tile_rotated = ([cosd(rotAngDeg) -sind(rotAngDeg) 0;
    sind(rotAngDeg) cosd(rotAngDeg) 0;
    0 0 1]*tile_pre_Rotate')';
tiles(1) = struct('nodes',tile_rotated(1:4,:),'nodeInds',[],'ndir',[]);
tiles(2) = struct('nodes',tile_rotated(5:8,:),'nodeInds',[],'ndir',[]);
tiles(3) = struct('nodes',tile_rotated(9:12,:),'nodeInds',[],'ndir',[]);
tiles(4) = struct('nodes',tile_rotated(13:16,:),'nodeInds',[],'ndir',[]);

% Doubly curved tiling
tiles(1) = struct('nodes',[0 0 0; 1 0 0; 1 1 0; 0 1 0],'nodeInds',[],'ndir',[]);
tiles(2) = struct('nodes',[1 1 0; 1+1/sqrt(2) 1 -1/sqrt(2); 1+1/sqrt(2) 2 -1/sqrt(2); 1 2 0],'nodeInds',[],'ndir',[]);
tiles(3) = struct('nodes',[0 2 0; 1 2 0; 1 2+1/sqrt(2) -1/sqrt(2); 0 2+1/sqrt(2) -1/sqrt(2)],'nodeInds',[],'ndir',[]);
tiles(4) = struct('nodes',[-1/sqrt(2) 1 -1/sqrt(2); 0 1 0; 0 2 0; -1/sqrt(2) 2 -1/sqrt(2)],'nodeInds',[],'ndir',[]);
%% Create connectivity graph
totNodesRaw = cat(1,tiles(:).nodes);        % array of all node coordinates in order of tiles before removing duplicates
lastTileNodesRawSep = zeros(length(tiles),1);   % array of the last index of each tile in order
% basenodes - coordinates of unique nodes
% ia - index of (one of) the original nodes that corresponds to that node in basenodes
% ic - index of basenodes that corresponds to that raw node
[basenodes, ia, ic] = unique(totNodesRaw,'rows','stable');
extrnodes = [];     % additional nodes added for rigid bodies and hinge constraints
isExtrBool = [];    % whether the nodes a rod is connected to are extra or base
tnodes = [];        % [b x 2] connectivity graph specifying the node indices at the end of the rod in unique basenode indices
for i = 1:length(tiles) % add necessary extra nodes and bars within each tile
    if rank(tiles(i).nodes-tiles(i).nodes(1,:)) ~= 2
        warning('Tile %i bad',i); % check coplanarity and non-collinearity
    end
    if i == 1
        lastTileNodesRawSep(i) = size(tiles(i).nodes,1);
    else
        lastTileNodesRawSep(i) = lastTileNodesRawSep(i-1) + size(tiles(i).nodes,1);
    end
    tiles(i).ndir = cross(tiles(i).nodes(3,:)-tiles(i).nodes(1,:),...
        tiles(i).nodes(2,:)-tiles(i).nodes(1,:));
    tiles(i).ndir = tiles(i).ndir*sign(tiles(i).ndir(3))/norm(tiles(i).ndir);
    [Lia,tiles(i).nodeInds] = ismember(tiles(i).nodes(:,:),basenodes,'rows'); % find the nodeInds for this tile's corners
    assert(all(Lia));
    tetraNodeInd = size(extrnodes,1)+1;     % Additional node to require rigidity
    extrnodes = [extrnodes; mean(tiles(i).nodes,1) + 0.1*tiles(i).ndir];
    tile_nodes = tiles(i).nodeInds;
    % Specify nodes and rods in base
    tnodes = [tnodes; tile_nodes(1) tile_nodes(2); tile_nodes(2) tile_nodes(3); tile_nodes(1) tile_nodes(3); ...
        tile_nodes(1) tetraNodeInd; tile_nodes(2) tetraNodeInd; tile_nodes(3) tetraNodeInd;]; % first tetra
    isExtrBool = [isExtrBool; false(3,2); repmat([false true],3,1)];
    for j = 4:length(tile_nodes)    % additional tetrahedrons
        tnodes = [tnodes; tile_nodes(1) tile_nodes(j); tile_nodes(j-1) tile_nodes(j);...
            tile_nodes(j) tetraNodeInd];
        isExtrBool = [isExtrBool; false(2,2); false true];
    end
end
%% Specify tile connections
% count - frequency of unique node (>1 means it's a joint)
% bin - index map from ic (raw nodes) to count
[count, ~, bin] = histcounts(ic,numel(ia));
assert(max(count)<=2);
% hinges = zeros(sum(count>1),3); % [tile1Ind tile2Ind nodeInd;]
hingeNodeInds = find(count>1);
for i = 1:sum(count>1)
    hingeNodeInd = hingeNodeInds(i);
    rawNodeInds = find(bin == hingeNodeInd);
    % want to take rawNodeInds and find the set of rawNodeInds in that respective tile
    tileInds = [find(lastTileNodesRawSep >= rawNodeInds(1),1) find(lastTileNodesRawSep >= rawNodeInds(2),1)]; % assumes only two tiles per hinge
    % Specify nodes and rods for hinges. Find nodes other than joint node
    t1Ind = tileInds(1);
    t2Ind = tileInds(2);
    t1Nodes = tiles(t1Ind).nodeInds;
    t1Nodes(t1Nodes == hingeNodeInd) = [];
    t2Nodes = tiles(t2Ind).nodeInds;
    t2Nodes(t2Nodes == hingeNodeInd) = [];
    t1ndir = tiles(t1Ind).ndir;
    t2ndir = tiles(t2Ind).ndir;
    if acos(t1ndir*t2ndir') < 1e5*eps % co-planar tiles
        extrnodes = [extrnodes; basenodes(hingeNodeInd,:) + 0.1*t1ndir];
        tnodes = [tnodes; hingeNodeInd length(extrnodes); t1Nodes(1) length(extrnodes); t1Nodes(2) length(extrnodes); ... % Link hinge joint to 2 other tile 1 points
            hingeNodeInd length(extrnodes); t2Nodes(1) length(extrnodes); t2Nodes(2) length(extrnodes)]; ... % Link hinge joint to 2 other tile 2 points
            isExtrBool = [isExtrBool; repmat([false true],3,1); ...
            repmat([false true],3,1)];
    else % not co-planar tiles
        extrnodes = [extrnodes; basenodes(hingeNodeInd,:) + 0.1*t1ndir; basenodes(hingeNodeInd,:) + 0.1*t2ndir];
        tnodes = [tnodes; hingeNodeInd length(extrnodes)-1; t1Nodes(1) length(extrnodes)-1; t1Nodes(2) length(extrnodes)-1; ... % Link hinge joint to 2 other tile 1 points
            hingeNodeInd length(extrnodes); t2Nodes(1) length(extrnodes); t2Nodes(2) length(extrnodes); ... % Link hinge joint to 2 other tile 2 points
            length(extrnodes) length(extrnodes)-1]; ... % link two external hinge joints
            isExtrBool = [isExtrBool; repmat([false true],3,1); ...
            repmat([false true],3,1); ...
            true true];
    end
end
nodes = [basenodes; extrnodes];
tnodes = tnodes + length(basenodes)*isExtrBool; % adjust so extranodes have greater index
for i = 1:length(tnodes) % re-arrange to make index increase
    if tnodes(i,1)>tnodes(i,2)
        temp = tnodes(i,1);
        tnodes(i,1) = tnodes(i,2);
        tnodes(i,2) = temp;
    end
end
%% Construct truss equilibrium matrix
% single hinge
% nodes = [0 0 0;...
%     1 0 0;...
%     1 2 0;...
%     0 2 0]; % [j x 3] coordinates of joints
% tnodes = [1 2;
%     2 3;
%     3 4;
%     1 4;
%     1 3]; % [b x 2] nodes each beam connects to, increasing index
% 3 chain
% nodes = [0 0 0;...
%     1 0 0;...
%     1 2 0;...
%     0 2 0
%     1 3 0]; % [j x 3] coordinates of joints
% tnodes = [1 2;
%     2 3;
%     3 4;
%     1 4;
%     1 3
%     3 5
%     4 5]; % [b x 2] nodes each beam connects to, increasing index
% tetrahedron
% nodes = [0 0 0;...
%     1 0 0;...
%     1 2 0;...
%     0.75 0.5 1]; % [j x 3] coordinates of joints
% tnodes = [1 2;
%     2 3;
%     1 3;
%     1 4;
%     2 4;
%     3 4]; % [b x 2] nodes each beam connects to, increasing index
% 2 tetrahedra
% nodes = [0 0 0;...
%     1 0 0;...
%     1 2 0;...
%     0.75 0.5 1
%     1 2 2.5
%     0.25 1.5 2]; % [j x 3] coordinates of joints
% tnodes = [1 2;
%     2 3;
%     1 3;
%     1 4;
%     2 4;
%     3 4;
%     3 5;
%     4 5;
%     3 6;
%     4 6;
%     5 6]; % [b x 2] nodes each beam connects to, increasing index
% 2 tetrahedra connected by 2R joint
% nodes = [0 0 0;...
%     1 0 0;...
%     1 2 0;...
%     0.75 0.5 1
%     0.9 0.7 1.5
%     1 2 2.5
%     0.25 1.5 2]; % [j x 3] coordinates of joints
% tnodes = [1 2;
%     2 3;
%     1 3;
%     1 4;
%     2 4;
%     3 4;
%     4 5;
%     3 5;
%     3 6;
%     3 7;
%     5 6;
%     5 7
%     6 7]; % [b x 2] nodes each beam connects to, increasing index

assert(nodes(1,1)==0  && nodes(1,2)==0 && nodes(1,3)==0); % first node origin
assert(nodes(2,1)~=0  && nodes(2,2)==0 && nodes(2,3)==0); % second node x axis
assert(nodes(3,2)~=0 && nodes(3,3)==0);  % third node xy plane
assert(all(tnodes(:,1)<tnodes(:,2)));
tdirs = nodes(tnodes(:,2),:) - nodes(tnodes(:,1),:); % [b x 3] unit directions of rods from low to high index
tdirs = tdirs./vecnorm(tdirs,2,2);
j = length(nodes);  % # of nodes
b = length(tdirs);  % # of rodes
A = zeros(3*j,b); % [3j x b], rows 1-3 correspond to node 1 x,y,z
for i = 1:b
    A(3*tnodes(i,1)-2:3*tnodes(i,1),i) = A(3*tnodes(i,1)-2:3*tnodes(i,1),i) + tdirs(i,:)';
    A(3*tnodes(i,2)-2:3*tnodes(i,2),i) = A(3*tnodes(i,2)-2:3*tnodes(i,2),i) - tdirs(i,:)';
end
%% SVD of equilibrium matrix
[U,V,W] = svd(A);
r = rank(A);
k = 6; % kinematic constraints for rigid motion
m = 3*j-k-r;  % mobility including infinitesimal mechanisms
s = b - r;    % states of selfstress
Umech = U(:,r+1:end); % inextensional mechanisms
Umech_rigid_trans = [repmat([1;0;0],j,1) repmat([0;1;0],j,1) repmat([0;0;1],j,1)];
Umech_rigid_rot = zeros(3*j,3);
for i = 1:j
    Umech_rigid_rot(3*i-2:3*i,1) = [0; -nodes(i,3); nodes(i,2)];
    Umech_rigid_rot(3*i-2:3*i,2) = [-nodes(i,3); 0; nodes(i,1)];
    Umech_rigid_rot(3*i-2:3*i,3) = [-nodes(i,2); nodes(i,1); 0];
end
Umech_rigid = [Umech_rigid_trans Umech_rigid_rot];
% Othogonalize rigid motion
Umech_rigid_ortho_norm = zeros(size(Umech_rigid));
for i1 = 1:6
    Umech_rigid_ortho_norm(:,i1) = Umech_rigid(:,i1);
    for i2 = 1:i1-1
        Umech_rigid_ortho_norm(:,i1) = Umech_rigid_ortho_norm(:,i1) ...
            - (Umech_rigid_ortho_norm(:,i1)'*Umech_rigid_ortho_norm(:,i2))*Umech_rigid_ortho_norm(:,i2);
    end
    Umech_rigid_ortho_norm(:,i1) = Umech_rigid_ortho_norm(:,i1)/norm(Umech_rigid_ortho_norm(:,i1));
end
% Check they are orthogonal
for i1 = 1:6
    for i2 = 1:6
        if i1 ~= i2
            assert(abs(Umech_rigid_ortho_norm(:,i1)'*Umech_rigid_ortho_norm(:,i2))<10*eps)
            %             if (abs(Umech_rigid_ortho_norm(:,i1)'*Umech_rigid_ortho_norm(:,i2))>10*eps)
            %                 'hi'
            %             end
        else
            assert(abs(Umech_rigid_ortho_norm(:,i1)'*Umech_rigid_ortho_norm(:,i2)-1)<10*eps)
        end
    end
end
% Orthogonalize internal mechanisms with rigid mechanisms and themselves
Umech_internal = zeros(size(Umech));    % internal mechanisms
for i1 = 1:size(Umech_internal,2)
    Umech_internal(:,i1) = Umech(:,i1);
    for i2 = 1:6
        Umech_internal(:,i1) = Umech_internal(:,i1) ...
            - (Umech(:,i1)'*Umech_rigid_ortho_norm(:,i2))*Umech_rigid_ortho_norm(:,i2); % shifts so resultant has zero sum x translation
    end
    for i2 = 1:i1-1
        Umech_internal(:,i1) = Umech_internal(:,i1) ...
            - (Umech(:,i1)'*Umech_internal(:,i2))*Umech_internal(:,i2);
    end
    if norm(Umech_internal(:,i1)) > 1e5*eps
        Umech_internal(:,i1) = Umech_internal(:,i1)/norm(Umech_internal(:,i1));
    else
        Umech_internal(:,i1) = zeros(3*j,1);
    end
end
Umech_internal = Umech_internal(:,vecnorm(Umech_internal)>10*eps);
if m ~= rank(Umech_internal)
    warning('Mobility not equal to SVD calculation. rank = %i. May have infinitesimal mechanism',rank(Umech_internal));
end
%% Plot visualization
figure
scatter3(nodes(:,1),nodes(:,2),nodes(:,3));
hold on
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal
% for i = 1:length(tnodes) % Plot load for self-stress
%     load = U(i,end);
%     if load > 0
%         plot3([nodes(tnodes(i,1),1) nodes(tnodes(i,2),1)]',...
%             [nodes(tnodes(i,1),2) nodes(tnodes(i,2),2)]',...
%             [nodes(tnodes(i,1),3) nodes(tnodes(i,2),3)]','k','LineWidth',load*10);
%     else
%         plot3([nodes(tnodes(i,1),1) nodes(tnodes(i,2),1)]',...
%             [nodes(tnodes(i,1),2) nodes(tnodes(i,2),2)]',...
%             [nodes(tnodes(i,1),3) nodes(tnodes(i,2),3)]','r','LineWidth',-load*10);
%     end
% end
plot3([nodes(tnodes(:,1),1) nodes(tnodes(:,2),1)]',...
    [nodes(tnodes(:,1),2) nodes(tnodes(:,2),2)]',...
    [nodes(tnodes(:,1),3) nodes(tnodes(:,2),3)]','k');
for i = 1:length(tiles)
    fill3(tiles(i).nodes(:,1), tiles(i).nodes(:,2), tiles(i).nodes(:,3),'blue','FaceAlpha',.3);
end
set(gca,'ColorOrderIndex',1)
for mech_ind = 1:rank(Umech_internal)
    mech_plot = Umech_internal(:,mech_ind);
    % Set node 1 displacement to 0
    mech_plot_based = mech_plot - Umech_rigid(:,1)*mech_plot(1);
    mech_plot_based = mech_plot_based - Umech_rigid(:,2)*mech_plot(2);
    mech_plot_based = mech_plot_based - Umech_rigid(:,3)*mech_plot(3);
    % Set node 2 displacement restricted to the line from the origin to it
    % rotate about y to get Umech_internal(6,:) to 0 by -mech_plot_based(6,:)/nodes(2,1)
    if nodes(2,1) == 0
        rotAng = 0;
    else
        rotAng = mech_plot_based(6)/nodes(2,1);
    end
    for i = 1:j
        mech_plot_based(3*i-2:3*i) = mech_plot_based(3*i-2:3*i) - [-nodes(i,3); 0; nodes(i,1)]*rotAng;
    end
    % rotate about z to get Umech_internal(5,:) to 0 by -mech_plot_based(5,:)/nodes(2,1)
    % get Umech_internal(3:5,:) in line with nodes(2,:)
    if nodes(2,1) == 0
        rotAng = 0;
    else
        rotAng = mech_plot_based(5)/nodes(2,1);
    end
    for i = 1:j
        mech_plot_based(3*i-2:3*i) = mech_plot_based(3*i-2:3*i) - [-nodes(i,2); nodes(i,1); 0]*rotAng;
    end
    % rotate about x to get Umech_internal(9,:) to 0 by -mech_plot_based(9,:)/nodes(3,2)
    if nodes(3,2) == 0
        rotAng = 0;
    else
        rotAng = mech_plot_based(9,:)/nodes(3,2);
    end
    for i = 1:j
        mech_plot_based(3*i-2:3*i) = mech_plot_based(3*i-2:3*i) - [0; -nodes(i,3); nodes(i,2)]*rotAng;
    end
%     quiver3(nodes(:,1), nodes(:,2), nodes(:,3), mech_plot_based(1:3:end), mech_plot_based(2:3:end), mech_plot_based(3:3:end));
    quiver3(nodes(1:size(basenodes,1),1), nodes(1:size(basenodes,1),2), nodes(1:size(basenodes,1),3),... % plot velocity only on base nodes
        mech_plot_based(1:3:3*size(basenodes,1)), mech_plot_based(2:3:3*size(basenodes,1)), mech_plot_based(3:3:3*size(basenodes,1)));
end
title(sprintf('Mobility = %i',m));