%% Make sections
DUMP_IMG = true;
OUTPUT_DIR = './img/';
% list of coronal slices
c_slices = unique(round(cs.y2i(E.xyz_entry(:,2))))';


h.fig_cslice = findobj('tag', 'fig_cslice');
if isempty(h.fig_cslice), h.fig_cslice = figure('color','w','tag', 'fig_cslice', 'Position', [100 100 800 700]); end
h.ax_cslice = findobj('tag', 'ax_cslice');
if isempty(h.ax_cslice), h.ax_cslice = axes('Parent', h.fig_cslice, 'tag', 'ax_cslice'); end


 
vol = V.phy;

for c = c_slices
    disp(c)
    isright = cs.x2i(E.xyz_entry(:,1)) <= cs.nx/2;
    id = round(cs.y2i(E.xyz0(:,2)))==c & isdeep & isright;
    is = round(cs.y2i(E.xyz0(:,2)))==c & ishallow & isright;
    
    imagesc(cs.xlim, cs.zlim, vol(:,:,c), 'parent', h.ax_cslice), colormap('bone');
    set(h.ax_cslice,'DataAspectRatio',[1 ,1,1], 'nextplot', 'add')
    
    lineplot = @(xyz0,xyz1,n) flatten([xyz0(:,n) xyz1(:,n) xyz1(:,n).*NaN]');
    col = get(gca,'colororder');
    plot((lineplot(E.xyz0(id,:), E.xyz_(id,:),1)),...
        (lineplot(E.xyz0(id,:), E.xyz_(id,:),3)), 'linewidth',2, 'color', col(1,:))
    plot((lineplot(E.xyz_entry(id,:), E.xyz_exit(id,:),1)),...
        (lineplot(E.xyz_entry(id,:), E.xyz_exit(id,:),3)),'-.', 'linewidth',1, 'color', col(1,:))
    plot((lineplot(E.xyz0(is,:), E.xyz_(is,:),1)),...
        (lineplot(E.xyz0(is,:), E.xyz_(is,:),3)), 'linewidth',2, 'color', col(2,:))
    plot((lineplot(E.xyz_entry(is,:), E.xyz_exit(is,:),1)),...
        (lineplot(E.xyz_entry(is,:), E.xyz_exit(is,:),3)),'-.', 'linewidth',1, 'color', col(2,:))
    set(h.ax_cslice, 'nextplot', 'replace')
    
    if DUMP_IMG
        if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end
        im_name = ['coronal_slice_' num2str(cs.i2y(c)*1e6, '%06.0fnm')];
        disp(im_name)
        if exist([OUTPUT_DIR im_name], 'file'), continue, end
        print(h.fig_cslice, [OUTPUT_DIR im_name], '-dpng')
        pause(0.05) % this is to stabilize as it sometimes fails
    end
end

disp('done')