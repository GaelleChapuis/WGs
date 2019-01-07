% 3D plot overlay
min_rec_length_mm = 2.5; % if probe insertion is smaller than than, discard
ylims = [-.009205 .004288]; % Antero Posterior selection (to remove OB and spine) WAXHOLM
mid_line_exclusion_mm = 0.4; % get away from mid-line vascular system

try delete(pl), end

% first take only the right side to have insertions towards medial plane
esel= cs.x2i(E.xyz_entry(:,1)) <= cs.nx/2;
esel = esel & E.rec_length > (min_rec_length_mm/1000); % prune according to the recording length
esel = esel & between(E.xyz_entry(:,2), ylims); % Remove olphactory bulb
% remove insertions too close to the midline vascular system
esel = esel & ~between(E.xyz_entry(:,1), cs.i2x(cs.nx/2)+[-1 1].*mid_line_exclusion_mm/1000);

T = [E.xyz_entry E.theta.*180/pi E.phi E.length.*1000 E.rec_length.*1000];
csv_head = ['X, Y, Z, Theta, Phi, Length_mm, Rec_length_mm' char(10)];
T = sortrows(T(esel,:), [2 1 4]);
CSV_IMPLANTS = 'implantations.csv';

fid = fopen(CSV_IMPLANTS,'w+'); fwrite(fid, csv_head); fclose(fid);
dlmwrite('implantations.csv', T,'-append')


% now display electrodes
col = get(gca,'colororder');
hold on, ie = E.theta==10/180*pi & esel;
pl(1) = plot3((E.xyz_entry(ie,1)), (E.xyz_entry(ie,2)), (E.xyz_entry(ie,3)), 'k*');          
pl(2) = plot3((flatten([E.xyz0(ie,1) E.xyz_(ie,1) E.xyz0(ie,1).*NaN ]')) ,...
              (flatten([E.xyz0(ie,2) E.xyz_(ie,2) E.xyz0(ie,2).*NaN ]')), ...
              (flatten([E.xyz0(ie,3) E.xyz_(ie,3) E.xyz0(ie,3).*NaN ]')), ...
               'color', col(4,:), 'parent', h.ax);
hold on, ie = E.theta==20/180*pi & esel;
pl(3) = plot3((flatten([E.xyz0(ie,1) E.xyz_(ie,1) E.xyz0(ie,1).*NaN ]')) ,...
              (flatten([E.xyz0(ie,2) E.xyz_(ie,2) E.xyz0(ie,2).*NaN ]')), ...
              (flatten([E.xyz0(ie,3) E.xyz_(ie,3) E.xyz0(ie,3).*NaN ]')),...
              'color', col(5,:), 'parent', h.ax);
set(pl,'linewidth',1.5)



