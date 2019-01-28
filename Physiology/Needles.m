function varargout = Needles(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Needles_OpeningFcn, ...
                   'gui_OutputFcn',  @Needles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function Needles_OpeningFcn(hobj, evt, h, varargin)
h.output = hobj;
guidata(hobj, h);

function varargout = Needles_OutputFcn(hobj, evt, h) 
logo = logo_ibl('square');
h.im_logo = image(h.axes_logo, logo);
axis(h.axes_logo, 'off');
% init preferences: when needed, to automatic check to add new prefs
default = struct('path_atlas', pwd,  'lateralize', true);
pfile = [fileparts(mfilename('fullpath')) filesep 'needles_param.json'];
if ~exist(pfile, 'file')
    io.write.json(pfile, default)
end
h.pref = io.read.json(pfile);
% wrap-up and save in handles
guidata(hobj, h);
varargout{1} = h.output;

function menu_file_loadatlas_Callback(hobj, evt, h)
[D.V, ~, D.cs, D.labels] = brainatlas.get_dsurqe(h.pref.path_atlas, 'lateralize', h.pref.lateralize, 'display', false);
lims.ap_lims = [-.009205 .004288]; % Antero Posterior selection (to remove OB and spine) WAXHOLM
lims = struct('ap_lims', [-0.005177 0.005503], 'ml_lims', [-0.004 0.004]); % for COG vertex, need to update once Bregma is chosen
% X: ML (pitch), 2nd_dim (left - right +)
% Y: AP (roll), 3d_dim +-
% Z: DV (yaw), 1st_dim -+
% VOL(Z, X, Y)
[D.E, D.S] = insert_electrodes(D.V, D.cs, lims);

% Create all the objects depending on the top axes
h.im_top = imagesc(D.cs.yscale, D.cs.xscale, D.S.top', 'Parent', h.axes_top);
set(h.axes_top, 'ydir', 'normal','DataAspectRatio',[1 1 1], 'NextPlot', 'add')
set(h.axes_top,'UserData', round(D.cs.y2i(0))) % WIndow motion callback
xlabel(h.axes_top, 'AP'), ylabel(h.axes_top, 'ML')
colormap(h.axes_top, 'bone')
h.pl_top_apline = plot(h.axes_top, [0 0], [D.cs.xlim], 'color', color_from_index(2));
set(h.pl_top_apline, 'ButtonDownFcn', @pl_top_apline_ButtonDownFcn)
h.pl_top_electrodes = plot(h.axes_top, D.E.xyz_entry(D.E.esel,2), D.E.xyz_entry(D.E.esel,1), '.', ...
'ButtonDownFcn', @pl_top_electrodes_ButtonDownFcn, 'MarkerSize', 10, 'Color', color_from_index(1));
h.pl_top_current_elec = plot(h.axes_top, NaN, NaN, '*', 'color', 'm', 'MarkerSize',12);
% Create all the objects depending on the contrast axes
h.im_phy = imagesc(D.cs.xscale, D.cs.zscale,  D.V.phy(:,:,round(D.cs.y2i(0))), 'Parent', h.axes_phy);
set(h.axes_phy, 'DataAspectRatio',[1 1 1], 'NextPlot', 'add')
h.pl_phy_xr = plot(h.axes_phy, [D.cs.xlim NaN 0 0], [0 0 NaN D.cs.ylim], 'Color', color_from_index(2));
h.pl_phy_zone = plot(h.axes_phy, NaN, 0,'.','MarkerSize',3,'Color',color_from_index(5));
h.pl_phy_electrodes(1) = plot(h.axes_phy, NaN, NaN, 'linewidth', 2, 'Color', color_from_index(3));
h.pl_phy_electrodes(2) = plot(h.axes_phy, NaN, NaN, 'linewidth', 2, 'Color', color_from_index(4));
h.pl_phy_electrodes_traj(1) = plot(h.axes_phy, NaN, NaN, '--', 'Color', color_from_index(3)); % Full electrode trajectory
h.pl_phy_electrodes_traj(2) = plot(h.axes_phy, NaN, NaN, '--', 'Color', color_from_index(4)); % Full electrode trajectory
h.pl_phy_current_elec = plot(h.axes_phy, NaN, NaN, '*', 'color', 'm', 'MarkerSize',12);
xlabel(h.axes_phy, 'ML'), ylabel(h.axes_phy, 'DV')
colormap(h.axes_phy, 'bone')
% Create all the objects depending on the label axes
h.im_lab = imagesc(D.cs.xscale, D.cs.zscale,  D.V.lab(:,:,round(D.cs.y2i(0))), 'Parent', h.axes_label);
set(h.axes_label, 'DataAspectRatio',[1 1 1], 'NextPlot', 'add')
h.pl_lab_xr = plot(h.axes_label, [D.cs.xlim NaN 0 0], [0 0 NaN D.cs.ylim], 'Color', color_from_index(2));
h.pl_lab_zone = plot(h.axes_label, NaN, 0,'.','MarkerSize',3,'Color',color_from_index(5));
h.pl_lab_electrodes(1) = plot(h.axes_label, NaN, NaN, 'linewidth', 2, 'Color', color_from_index(3));
h.pl_lab_electrodes(2) = plot(h.axes_label, NaN, NaN, 'linewidth', 2, 'Color', color_from_index(4));
h.pl_lab_electrodes_traj(1) = plot(h.axes_label, NaN, NaN, '--', 'Color', color_from_index(3)); % Full electrode trajectory
h.pl_lab_electrodes_traj(2) = plot(h.axes_label, NaN, NaN, '--', 'Color', color_from_index(4)); % Full electrode trajectory
h.pl_lab_current_elec = plot(h.axes_label, NaN, NaN, '*', 'color', 'm', 'MarkerSize',12);
xlabel(h.axes_label, 'ML'), ylabel(h.axes_label, 'DV')
colormap(h.axes_label, 'bone')
% Init graphic states for callbacks
set([h.pl_phy_xr, h.pl_lab_xr, h.pl_phy_zone, h.pl_lab_zone], 'Visible', 'Off')
set(h.fig_main,'WindowButtonMotionFcn', {@fig_main_WindowButtonMotionFcn, h})
% prevents from re-loading the Atlas for the time being
set(h.menu_file_loadatlas, 'enable', 'off')
guidata(h.fig_main, h)
setappdata(h.fig_main, 'Data', D)


function pl_top_electrodes_ButtonDownFcn(hobj, evt)
h = guidata(hobj);
D = getappdata(h.fig_main, 'Data');
ie = find((abs(evt.Source.XData-evt.IntersectionPoint(1)) < 1e-9) & ...
(abs(evt.Source.YData-evt.IntersectionPoint(2)) < 1e-9));
% Find the electrode index from the plot to the data structure
esel = find(D.E.esel);
ie = esel(ie);
ap_current =D.E.xyz0(ie(1),2);
Update_Slices(h.fig_main, [], ap_current);
set(h.axes_top,'UserData', round(D.cs.y2i(ap_current(1))))
Update_txt_electrodes(hobj, ie);
set([ h.pl_lab_current_elec, h.pl_phy_current_elec],'Visible', 'on',...
    'xdata', D.E.xyz_entry(ie,1), 'ydata', D.E.xyz_entry(ie,3))
set( h.pl_top_current_elec,'Visible', 'on',...
    'xdata', D.E.xyz_entry(ie,2), 'ydata', D.E.xyz_entry(ie,1))


function pl_top_apline_ButtonDownFcn(hobj, evt)
h = guidata(hobj);
set(h.axes_top,'UserData', NaN)
set(h.fig_main, 'WindowButtonUpFcn', @pl_top_apline_ButtonUpFcn)
set([h.pl_lab_zone h.pl_phy_zone], 'xdata', NaN, 'ydata', NaN)
set(h.txt_labels, 'String', '')

function pl_top_apline_ButtonUpFcn(hobj, evt)
h = guidata(hobj);
D = getappdata(hobj, 'Data');
ap_current = get(h.pl_top_apline, 'xdata');
Update_Slices(h.fig_main, [], ap_current)
set(h.axes_top,'UserData', round(D.cs.y2i(ap_current(1))))
set(h.fig_main, 'WindowButtonUpFcn', '')

function fig_main_WindowButtonMotionFcn(hobj, evt, h)
% Next gets executed if over slice or label axes
xy = get(h.axes_phy, 'CurrentPoint');
if between(xy(1), get(h.axes_phy,'xlim')) && between(xy(3), get(h.axes_phy,'ylim'))
    Update_CrossHairs(hobj, evt, xy([1 3]))
    Update_Labels(hobj, evt, xy([1 3]))
    return
end
xy = get(h.axes_label, 'CurrentPoint');
if between(xy(1), get(h.axes_label,'xlim')) && between(xy(3), get(h.axes_label,'ylim'))
    Update_CrossHairs(hobj, evt, xy([1 3]))
    Update_Labels(hobj, evt, xy([1 3]))
    return
end
% stuff gets executed if over top axes
xy = get(h.axes_top, 'CurrentPoint');
if between(xy(1), get(h.axes_top,'xlim')) && between(xy(3), get(h.axes_top,'ylim'))
    set([h.pl_phy_xr, h.pl_lab_xr], 'Visible', 'Off')
    Update_txt_xyz(hobj, xy(2), xy(1), NaN)
    if isnan(get(h.axes_top, 'UserData'))
        Update_Slices(h.fig_main, [], xy(1));
    end
    return
end
set([h.pl_phy_xr, h.pl_lab_xr], 'Visible', 'Off')
set(h.txt_xyz,'Visible', 'off');


function Update_CrossHairs(hobj, evt, ml_dv)
h = guidata(hobj);
D = getappdata(hobj, 'Data');
set([h.pl_phy_xr, h.pl_lab_xr, h.pl_phy_zone, h.pl_lab_zone], 'Visible', 'On')
set( [h.pl_lab_xr, h.pl_phy_xr], 'xdata', [D.cs.xlim NaN ml_dv([1 1])],...
    'ydata', [ml_dv([2 2]) NaN D.cs.ylim]);
ap = D.cs.i2y(get(h.axes_top, 'UserData'));
Update_txt_xyz(h.txt_xyz, ml_dv(1), ap, ml_dv(2));


function Update_txt_electrodes(hobj, ie)
h = guidata(hobj);
D = getappdata(h.fig_main, 'Data');
set(h.txt_electrode, 'Visible', 'on', 'String', { 'ELECTRODE:', '', ...
    ['ML Line: ' num2str(D.E.Line(ie(1)), '%04.0f')], ...
    ['AP Point: ' num2str(D.E.Point(ie(1)), '%04.0f')], ...
    ['AP: ' num2str(D.E.xyz_entry(ie(1),2)*1e6, '%06.0f (um)')],...
    ['ML: ' num2str(D.E.xyz_entry(ie(1),1)*1e6, '%06.0f (um)')]});
%     ['Depths (shallow): ' num2str(D.E.length(ie(1)',1)*1e6, '%06.0f (um)')],...
%     ['Depths (deep): ' num2str(D.E.length(ie(2)',1)*1e6, '%06.0f (um)')]},...
    

function Update_txt_xyz(hobj, ml, ap, dv)
h = guidata(hobj);
set(h.txt_xyz, 'String', {['ML:', num2str(ml.*1e6, '%06.0f (um)')],'',...
                          ['AP:', num2str(ap.*1e6, '%06.0f (um)')],'',...
                          ['DV:', num2str(dv.*1e6, '%06.0f (um)')]},...
                          'Visible', 'on')


function Update_Labels(hobj, evt, ml_dv)
h = guidata(hobj);
D = getappdata(hobj, 'Data');
yi = get(h.axes_top,'UserData');
lab = D.V.lab(  round(D.cs.z2i(ml_dv(2))), round(D.cs.x2i(ml_dv(1))), yi);
if lab==0, return, end
[x,z]= find(D.V.lab(:,:,yi)== lab);
set([h.pl_phy_zone h.pl_lab_zone], 'xdata', D.cs.i2x(z), 'ydata', D.cs.i2z(x), 'visible', 'on')
structure = D.labels.name{D.labels.index==lab};
try
set(h.txt_labels, 'String', structure)
catch
    disp(structure)
end


function Update_Slices(hobj, evt, ap)
D = getappdata(hobj, 'Data');
h = guidata(hobj);
set(h.pl_top_apline, 'Xdata', ap([1 1]));
ap_slice = round(D.cs.y2i( ap(1)) );
if between(ap_slice,[1 D.cs.ny])
    set(h.im_phy, 'CData', D.V.phy(:,:,ap_slice))
    set(h.im_lab, 'CData', D.V.lab(:,:,ap_slice))
end
set([h.pl_lab_zone, h.pl_phy_zone], 'visible', 'off')
set([h.pl_top_current_elec h.pl_phy_current_elec h.pl_lab_current_elec], 'visible', 'off')
% Find the electrodes from the closest coronal plane
[d, ie] = min(abs(ap(1) -  D.E.xyz0(:,2)));
ie = D.E.xyz0(:,2)==D.E.xyz0(ie,2);
% Plot Electrodes
lineplot = @(xyz0,xyz1,n) flatten([xyz0(:,n) xyz1(:,n) xyz1(:,n).*NaN]');
% plot 10 degres insertions, active shank and full track
i1 = D.E.esel & ie & D.E.theta == 10*pi/180;
set([h.pl_phy_electrodes(1) h.pl_lab_electrodes(1)],...
    'xdata', lineplot(D.E.xyz0(i1,:), D.E.xyz_(i1,:),1),...
    'ydata', lineplot(D.E.xyz0(i1,:), D.E.xyz_(i1,:),3))
set([h.pl_phy_electrodes_traj(1) h.pl_lab_electrodes_traj(1)],...
    'xdata', lineplot(D.E.xyz_entry(i1,:), D.E.xyz_exit(i1,:),1),...
    'ydata', lineplot(D.E.xyz_entry(i1,:), D.E.xyz_exit(i1,:),3))
% plot 20 degres insertions, active shank and full track
i2 = D.E.esel & ie & D.E.theta == 20*pi/180;
set([h.pl_phy_electrodes(2) h.pl_lab_electrodes(2)],...
    'xdata', lineplot(D.E.xyz0(i2,:), D.E.xyz_(i2,:),1),...
    'ydata', lineplot(D.E.xyz0(i2,:), D.E.xyz_(i2,:),3))
set([h.pl_phy_electrodes_traj(2) h.pl_lab_electrodes_traj(2)],...
    'xdata', lineplot(D.E.xyz_entry(i2,:), D.E.xyz_exit(i2,:),1),...
    'ydata', lineplot(D.E.xyz_entry(i2,:), D.E.xyz_exit(i2,:),3))
