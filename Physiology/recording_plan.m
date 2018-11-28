%% All units SI
% http://help.brain-map.org/display/mousebrain/API#API-DownloadImages
% http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/ara_nissl/
file_label_csv = '/s0/BrainAtlas/structure_tree_safe_2017.csv';
res = 50;
switch res
    case 50
        nrd_file_annotations =  '/s0/BrainAtlas/annotation_50.nrrd';
        nrd_file_nissl =  '/s0/BrainAtlas/ara_nissl_50.nrrd';
    case 100
        nrd_file_annotations =  '/s0/BrainAtlas/annotation_100.nrrd';
        nrd_file_nissl =  '/s0/BrainAtlas/ara_nissl_100.nrrd';
end

[dx,dy,dz] = deal(res.*1e-6); % sampling interval of model (m)
labels = readtable(file_label_csv);
[av, hn] = io.read.nrrd(nrd_file_annotations );
[gs, hn] = io.read.nrrd(nrd_file_nissl );
% Organization in memory should reflect most frequent usage: coronal first, sagittal second and transverse last
% This is not important if the file is small, for a large file and memory mapping this will make the difference between usable and unusable
av = permute(av,[1,3,2]); 
gs = permute(gs,[1,3,2]); 

[id, ind] = unique(av(:));
V= av;
for m = 1:length(ind)
    V(av==id(m)) = ind(m)-1;
end

% X: ML (pitch), 2nd_dim
% Y: AP (roll), 3d_dim
% Z: DV (yaw), 1st_dim
% cerebral aqueduct
% ventricular system

% length of each dimension
nx =size(V,2); ny = size(V,3); nz = size(V,1);
% distance to indices
x2i = @(x) x./dx+1; y2i = @(y) y./dy+1; z2i = @(z) z./dz+1;
% indices to distance
i2x = @(m) (m-1).*dx; i2y = @(m) (m-1).*dy; i2z = @(m) (m-1).*dz; %
% relative position (ratio) to indices
rx2i = @(x) x.*(nx-1) +1; ry2i = @(y) y.*(ny-1)+1; rz2i = @(z) z.*(nz-1)+1; 
% relative position (ratio) to distance
rx2x = @(x) i2x(rx2i(x)); ry2y = @(x) i2y(ry2i(x)); rz2z = @(x) i2z(rz2i(x));

%%
close all
fv = isosurface(permute(V~=0,[3, 2, 1]),0.5);
h.fig = figure('Color','w'); h.p = patch(fv); h.ax = gca;
set(h.ax, 'DataAspectRatio',[1 1 1], 'zdir', 'reverse')
xlabel(h.ax, 'x'), ylabel(h.ax, 'y'), zlabel(h.ax, 'z')
h.p.FaceColor = 'red';
h.p.EdgeColor = 'none';
h.p.FaceAlpha = 0.7;
view(69,42);
camlight;

% extract surfaces of the top and bottom of the brain
[S.top, S.bottom] = deal(zeros(size(V,3), size(V,2)));
for s=1:size(V,3)
    [i1,i2] = find(diff(V(:,:,s)==0,1,1));
    S.top(s,:) = accumarray(i2, i1, [size(V,2) 1], @min, NaN);
    S.bottom(s,:) = accumarray(i2, i1, [size(V,2) 1], @max, NaN);
end
% figure
%  hold on, surf(S.top, 'EdgeColor', 'none'), 
%  hold on, surf(S.bottom, 'EdgeColor', 'none'), 

%% Electrodes dans le plan coronal
len_active_electrode = 3.8 *1e-3;
% len_active_electrode = 10 *1e-3;
len_electrode_tip = 1*1e-3;
len_electrode = 11*1e-3; % fixed electrode length(tip+shank)
phi = 0; % always in coronal plan
theta = 10/180*pi; % insertion is 10 degrees here
delec = 0.5*1e-3; % we start with a grid of half mm
% create a mesh of electrodes, extract the z coordinate from the brain surface
E = [];
[y_, x_] = meshgrid([0:delec:i2x(size(V,3))], [0:delec:i2x(size(V,2))]);
z_ = i2z(interp2(S.top, x2i(x_), y2i(y_)));
E.xyz_entry = [x_(~isnan(z_(:))), y_(~isnan(z_(:))), z_(~isnan(z_(:)))];
E.phi = E.xyz_entry(:,1).*0 + phi; % always in coronal plan
E.theta = E.phi.*0;
% every insertion site is used twice: 10 and 20 degrees
ne = size(E.xyz_entry,1);
E = structfun(@(x) repmat(x,2,1), E, 'UniformOutput', false);
E.theta(1:ne) = 10/180*pi; % shallow electrodes 10 degrees
E.theta(ne+1:end) = 20/180*pi; % deep electrodes 20 degrees
ne = size(E.xyz_entry,1);

 
% Find the electrodes path using polar coordinates
r = len_electrode; % fixed electrode length (tip+sites)
E.xyz_ = probe_sph2cart(r, E.theta, E.phi, E.xyz_entry);

% get a uniform sampling of points along the electrode paths to compute exit point
nr = ceil(len_electrode*4/dx); % number of sample points along the path
X_ = x2i((E.xyz_(:,1) - E.xyz_entry(:,1))*linspace(0,1,nr) + E.xyz_entry(:,1)) ;
Y_ = y2i((E.xyz_(:,2) - E.xyz_entry(:,2))*linspace(0,1,nr) + E.xyz_entry(:,2));
Z_ = z2i((E.xyz_(:,3) - E.xyz_entry(:,3))*linspace(0,1,nr) + E.xyz_entry(:,3));
E.v = interp3(V, X_, Z_, Y_, 'nearest'); %extract the indices from the volume
% if all points are 0 this is a dud but this shouldn't happen as each
% electrode has a point of entry on the top surface
assert(~any(all(E.v==0,2)))
[aa,bb] = find(diff(E.v==0,1,2)==1);
E.length = E.phi.*0;
E.length(aa) = bb./(nr-1).*r;
E.xyz_exit = probe_sph2cart(E.length, E.theta, E.phi, E.xyz_entry);

% compute active path part of electrode only
E.xyz0 = E.xyz_entry.*0;
% for shallow electrodes entrypoint + active length
ishallow = (E.theta == 10/180*pi);
E.xyz0(ishallow,:) = E.xyz_entry(ishallow,:);
E.rec_length = min(len_active_electrode, max(0, E.length - len_electrode_tip));
E.xyz_(ishallow,:) = probe_sph2cart(E.rec_length(ishallow), E.theta(ishallow), E.phi(ishallow), E.xyz0(ishallow,:));
% for deep electrodes start from bottom - tip and up to the rec length
isdeep = (E.theta == 20/180*pi);
E.xyz_(isdeep,:) = probe_sph2cart(-len_electrode_tip, E.theta(isdeep), E.phi(isdeep), E.xyz_exit(isdeep,:));
E.rec_length(isdeep) = min(len_active_electrode, max(0, E.length(isdeep) - len_electrode_tip));
E.xyz0(isdeep,:) = probe_sph2cart(-E.rec_length(isdeep), E.theta(isdeep), E.phi(isdeep), E.xyz_(isdeep,:)); 

%% 3D plot overlay
try delete(pl), end
% figure,
% hold on, surf(S.bottom, 'EdgeColor', 'None'); set(gca, 'DataAspectRatio',[1 1 1])
% hold on, surf(S.top, 'EdgeColor', 'None'); set(gca, 'DataAspectRatio',[1 1 1])
isleft = x2i(E.xyz_entry(:,1)) >= nx/2;
col = get(gca,'colororder');
hold on, ie = E.theta==10/180*pi & isleft ;
pl(1) = plot3(x2i(E.xyz_entry(ie,1)), y2i(E.xyz_entry(ie,2)), z2i(E.xyz_entry(ie,3)), 'k*');          
pl(2) = plot3(x2i(flatten([E.xyz0(ie,1) E.xyz_(ie,1) E.xyz0(ie,1).*NaN ]')) ,...
              y2i(flatten([E.xyz0(ie,2) E.xyz_(ie,2) E.xyz0(ie,2).*NaN ]')), ...
              z2i(flatten([E.xyz0(ie,3) E.xyz_(ie,3) E.xyz0(ie,3).*NaN ]')), 'color', col(4,:));
hold on, ie = E.theta==20/180*pi & isleft ;
pl(3) = plot3(x2i(flatten([E.xyz0(ie,1) E.xyz_(ie,1) E.xyz0(ie,1).*NaN ]')) ,...
              y2i(flatten([E.xyz0(ie,2) E.xyz_(ie,2) E.xyz0(ie,2).*NaN ]')), ...
              z2i(flatten([E.xyz0(ie,3) E.xyz_(ie,3) E.xyz0(ie,3).*NaN ]')), 'color', col(5,:));
set(pl,'linewidth',1.5)
          
%% Make sections
for c = unique(round(y2i(E.xyz_entry(:,2))))'
    disp(c)
isright = x2i(E.xyz_entry(:,1)) <= nx/2;
id = round(y2i(E.xyz0(:,2)))==c & isdeep & isright;
is = round(y2i(E.xyz0(:,2)))==c & ishallow & isright;

f = figure('color','w')
imagesc(gs(:,:,c)), colormap('bone'); 
set(gca,'DataAspectRatio',[1 ,1,1])
hold on,

lineplot = @(xyz0,xyz1,n) flatten([xyz0(:,n) xyz1(:,n) xyz1(:,n).*NaN]');
col = get(gca,'colororder');
plot(x2i(lineplot(E.xyz0(id,:), E.xyz_(id,:),1)),...
     z2i(lineplot(E.xyz0(id,:), E.xyz_(id,:),3)), 'linewidth',2, 'color', col(1,:))
plot(x2i(lineplot(E.xyz_entry(id,:), E.xyz_exit(id,:),1)),...
     z2i(lineplot(E.xyz_entry(id,:), E.xyz_exit(id,:),3)),'-.', 'linewidth',1, 'color', col(1,:))
plot(x2i(lineplot(E.xyz0(is,:), E.xyz_(is,:),1)),...
     z2i(lineplot(E.xyz0(is,:), E.xyz_(is,:),3)), 'linewidth',2, 'color', col(2,:))
plot(x2i(lineplot(E.xyz_entry(is,:), E.xyz_exit(is,:),1)),...
     z2i(lineplot(E.xyz_entry(is,:), E.xyz_exit(is,:),3)),'-.', 'linewidth',1, 'color', col(2,:))
 
 im_name = ['coronal_slice_' num2str(c,'%03.0f')];
 print(f, ['./img/' im_name], '-dpng')
end





%% Histogram


%% Compute the coverage map
% #NB from neuropixel paper, it is hard to see anuthing beyond 60um away
% from the probe. However we'll assume that 0.5mm is what we're looking for
% so linear score 1/2 @ 250um
fcn_coverage_score = @(r) max(0,1-r./(500.*1e-6)); % This is trivial but we may want to change it later without re-writing the code
%fcn_coverage_score([0,250,500,800].*1e-6) = [1, 0.5, 0, 0]

% now need to compute distance to electrode
c_slices = unique(round(y2i(E.xyz0(:,2))));





