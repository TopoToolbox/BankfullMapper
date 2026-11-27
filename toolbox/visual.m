function [final_Zlim,final_Zvar,final_Qlim,final_Qvar,residual_Z,residual_Q,res_Z,res_Q]=visual(HAR,DEM,SW,step,zlim,zvar,qlim,qvar,varargin)

%     map visualization of the bankfull geometry and discharge 
%
% Syntax
%
%    [final_Zlim,final_Zvar,final_Qlim,final_Qvar,residual_Z,residual_Q,res_Z,res_Q]=visual(HAR,DEM,SW,step,zlim,zvar,qlim,qvar)
%
% Description
%
%     This function allows to map the bankfull geometry and discharge based
%     on the results obtained by DETECT_PEAK and MANNINGEQ functions. In
%     detail, both peak elevation values and most probable values for each
%     section are used for the bankfull geometry visualization. Regarding
%     the discharge visualization, both the discharge computed with the
%     peak elevation values and that with the most probable elevation
%     values are used. The bankfull elevation and discharge residuals are
%     also computed by subtracting the most probable values with the peak values.
%
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     HAR       high resolution GRIDobj
%     SW        SWATHobj of the river profile
%     step      stepping in meters of the profiles
%     zlim      peak elevation above thalveg values (see DETECT_PEAK results)
%     zvar      most probable elevation values above thalveg values (see DETECT_PEAK results)
%     qlim      discharge values (computed with peak elevation values zlim using MANNINGEQ function)
%     qvar      discharge values (computed with the most probable elevation values zvar using MANNINGEQ function)
%
% Optional Input arguments
%
%     connected true (default): only the flooded cells connected
%               (8-neighbourhood) to the thalweg are retained, so that
%               isolated depressions within the swath (ditches, ponds,
%               oxbows) are not mapped. false: original behaviour
%     thalwegwindow half-width in meters of the window around the river
%               course where the thalweg of each transverse line of the
%               swath is searched as the lowest point (default 10; use the
%               same value used in BANKFULL, DETECT_PEAK and MANNINGSEQ)
%     nearest   true (default): each cell within half the swath width from
%               the river course is assigned to the nearest point of the
%               course, and then to the nearest profile along the course
%               (robust where transverse profiles overlap, e.g. in tight
%               meanders). false: original approach (along-course distance
%               interpolated from the SWATHobj with SWATHobj2GRIDobj)
%     plot      plotting the results
%
%     Each profile is applied to a reach centred on the profile itself
%     (boundaries halfway between consecutive profiles).
%
% Output arguments
% 
%     final_Zlim        GRIDobj of the REM with the values below the most probable peak elevation 
%     final_Zvar        GRIDobj of the REM with the variable z values 
%     final_Qlim        GRIDobj of the REM with discharge values with constant z 
%     final_Qvar        GRIDobj of the REM with discharge values with variable z
%     residual_Z        GRIDobj of the residual between the residual error between
%                       the peak elevation of each section and the most probable 
%                       elevation 
%     residual_Q        GRIDobj of the residual between the residual error between
%                       the discharge of each section extracted with constant z and the most probable 
%                       discharge 
%     res_Z             array of the elevation residual values    
%     res_Q             array of the discharge values extracted with constant z     
%
%
% 
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024, - updated 27. November, 2026
%
% Parse input
p = inputParser;
p.FunctionName = 'visual';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'HAR',@(x) isa(x,'GRIDobj'));
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParameter(p,'step',@(x) isnumeric(x))
addParameter(p,'zlim',@(x) isnumeric(x))
addParameter(p,'zvar',@(x) isnumeric(x))
addParameter(p,'qlim',@(x) isnumeric(x))
addParameter(p,'qvar',@(x) isnumeric(x))
addParameter(p,'connected',true,@(x) isscalar(x))
addParameter(p,'nearest',true,@(x) isscalar(x))
addParameter(p,'thalwegwindow',10,@(x) isnumeric(x) && isscalar(x) && x>=0)
addParameter(p,'plot',true);

parse(p,DEM,HAR,SW,varargin{:});

plt=p.Results.plot;
connected=logical(p.Results.connected);
nearest=logical(p.Results.nearest);
twin=p.Results.thalwegwindow;

%% SECTIONING OF STREAM PROFILE
dist_along=SW.distx;
dx = SW.dx;
dix = round(step/dx);
ix = 1:dix:length(dist_along);

% along-course distance of the profiles: each profile is applied to a
% reach centred on it, with boundaries halfway between consecutive
% profiles (class i of DX = section i)
sd = dist_along(ix);
sd = sd(:)';

% seeds for the connectivity with the channel: thalweg of each transverse
% line of the swath, i.e. its lowest point within thalwegwindow meters from
% the river course (the same rule used in BANKFULL, DETECT_PEAK, MANNINGSEQ)
Zw = SW.Z;
Zw(abs(SW.disty) > twin,:) = NaN;
Zw(isnan(Zw)) = Inf;
[zmin,rmin] = min(Zw,[],1);
cols = find(isfinite(zmin));
iseed = sub2ind(size(SW.Z),rmin(cols),cols);
TRC = HAR;
TRC.Z = false(size(HAR.Z));
itr = coord2ind(HAR,SW.X(iseed),SW.Y(iseed));
TRC.Z(itr(~isnan(itr))) = true;

if nearest
    % each cell within half the swath width from the river course is
    % assigned to the nearest point of the course (robust where the
    % transverse profiles overlap, e.g. in tight meanders)
    [DX,T] = nearestsection(HAR,SW,sd);
else
    % original approach: along-course distance interpolated from the swath
    DX = SWATHobj2GRIDobj(SW,HAR,'distx');
    if numel(sd) > 1
        edges = (sd(1:end-1)+sd(2:end))/2;
        DX = reclassify(DX,'definedintervals',edges(:));
    else
        DX.Z(~isnan(DX.Z)) = 1;
    end
    T=SWATHobj2GRIDobj(SW,HAR);
    T=T~=0;
end

% clip HAR in the SW area
REMc=clip(HAR,T);

%% reclassification and plotting of final_Zlim, final_Zvar, residual_Z and res_Z

%
final_Zlim=REMc;

for i =1:length(zlim)
    final_Zlim.Z(DX.Z==i)=final_Zlim.Z(DX.Z==i)<zlim(i);
end

M=final_Zlim==1;
if connected
    M.Z = imreconstruct(M.Z & TRC.Z, M.Z, 8);
end

for i =1:length(zlim)
    final_Zlim.Z(DX.Z==i)=zlim(i);
end

final_Zlim=clip(final_Zlim,M);

%
final_Zvar=REMc;

for i =1:length(zvar)
    final_Zvar.Z(DX.Z==i)=final_Zvar.Z(DX.Z==i)<zvar(i);
end

M=final_Zvar==1;
if connected
    M.Z = imreconstruct(M.Z & TRC.Z, M.Z, 8);
end

for i =1:length(zvar)
    final_Zvar.Z(DX.Z==i)=zvar(i);
end

final_Zvar=clip(final_Zvar,M);

% 
res_Z=zlim-zvar;

residual_Z=REMc;


for i =1:length(res_Z)
    residual_Z.Z(DX.Z==i)=res_Z(i);
    
end

M=~isnan(final_Zvar);
residual_Z=clip(residual_Z,M);

%% reclassification and plotting of final_Qlim, final_Qvar, residual_Q and res_Q

%
final_Qlim=final_Zlim;

for i =1:length(qlim)
    final_Qlim.Z(DX.Z==i)=qlim(i);
    
end

M=~isnan(final_Zlim);
final_Qlim=clip(final_Qlim,M);

%
final_Qvar=final_Zvar;

for i =1:length(qvar)
    final_Qvar.Z(DX.Z==i)=qvar(i);
    
end

M=~isnan(final_Zvar);
final_Qvar=clip(final_Qvar,M);

% 
res_Q=qlim-qvar;

residual_Q=REMc;


for i =1:length(res_Q)
    residual_Q.Z(DX.Z==i)=res_Q(i);
    
end

M=~isnan(final_Qvar);
residual_Q=clip(residual_Q,M);



%% plotting 
if plt

% plot final_Z
figure
imageschs(DEM,final_Zlim,'caxis',[0 max(final_Zlim)],'colorbarylabel','Elevation above thalveg (m)')

title('Final Zlim')

% plot final_Zvar
figure
imageschs(DEM,final_Zvar,'caxis',[0 max(final_Zvar)],'colorbarylabel','Elevation above thalveg (m)')

title('Final Zvar')

% plot final_Qlim
figure
imageschs(DEM,final_Qlim,'caxis',[0 max(final_Qlim)],'colorbarylabel','Discharge (m^{3} s^{-1})')

title('Final Q')

% plot final_Qvar
figure
imageschs(DEM,final_Qvar,'caxis',[0 max(final_Qvar)],'colorbarylabel','Discharge (m^{3} s^{-1})')

title('Final Q with variable z')

% plot residual_Z
figure
imageschs(DEM,residual_Z,'caxis',[-max(residual_Z) max(residual_Z)],'colormap','jet','colorbarylabel','Residual (m)')

title('Residual Z')   

% plot residual_Q
figure
imageschs(DEM,residual_Q,'caxis',[min(residual_Q) max(residual_Q)],'colormap','jet','colorbarylabel','Residual discharge (m^{3} s^{-1})')

title('Residual discharge') 


% plot res_Z and res_Q

figure
subplot(2,1,1)

plot(zlim,'k-')
hold on
plot(zvar,'r-')

ylabel('Peak elevation (m)')
xlabel('Section')

yyaxis right
plot(res_Z,'g')
yline(0)

ylabel('Residual (m)')

legend('Peak elevation','Most probable peak elevation','Residual')
  
subplot(2,1,2)
plot(qlim,'k-')
hold on
plot(qvar,'r-')

ylabel('Discharge (m^{3} s^{-1})')
xlabel('Section')

yyaxis right
plot(res_Q,'g')
yline(0)

ylabel('Residual discharge (m^{3} s^{-1})')

legend('Peak discharge','Most probable peak discharge','Residual')

end

end

%% ------------------------------------------------------------------------
function [DX,T] = nearestsection(HAR,SW,sd)
% section index of each cell (0 outside the swath) and swath footprint,
% assigning each cell to the nearest point of the river course

[gx,gy] = getcoordinates(HAR);
gx = gx(:); gy = gy(:);
dxg = gx(2)-gx(1);
dyg = gy(2)-gy(1);
hw  = SW.width/2;

% subgrid covering the swath
sx = [SW.X(:); SW.xy(:,1)];
sy = [SW.Y(:); SW.xy(:,2)];
v  = ~isnan(sx) & ~isnan(sy);
c  = round((sx(v)-gx(1))/dxg)+1;
r  = round((sy(v)-gy(1))/dyg)+1;
c1 = max(1,min(c)-1); c2 = min(numel(gx),max(c)+1);
r1 = max(1,min(r)-1); r2 = min(numel(gy),max(r)+1);
nr = r2-r1+1; nc = c2-c1+1;

% cells of the river course (to delimit the swath)
tr = round((SW.xy(:,2)-gy(1))/dyg)+1 - r1 + 1;
tc = round((SW.xy(:,1)-gx(1))/dxg)+1 - c1 + 1;
in = tr>=1 & tr<=nr & tc>=1 & tc<=nc;
if ~any(in)
    error('visual:trace','The river course does not overlap the HAR grid.')
end
TR = false(nr,nc);
TR(sub2ind([nr nc],tr(in),tc(in))) = true;

% swath footprint: cells within half the swath width from the river course
% (distance only: BWDIST does not return indices for images larger than
% 2^24 pixels)
D = bwdist(TR);
foot = double(D)*HAR.cellsize <= hw & ~isnan(HAR.Z(r1:r2,c1:c2));
clear D

% along-course distance of the nearest point of the river course for each
% cell of the footprint (processed in blocks to limit memory)
[fr,fc] = find(foot);
along = nan(numel(fr),1);
NS = KDTreeSearcher(SW.xy);
chunk = 1e6;
for k = 1:chunk:numel(fr)
    ii = k:min(k+chunk-1,numel(fr));
    idx = knnsearch(NS,[gx(c1+fc(ii)-1) gy(r1+fr(ii)-1)]);
    along(ii) = SW.distx(idx);
end

% nearest profile along the course
if numel(sd) > 1
    edges = [-inf, (sd(1:end-1)+sd(2:end))/2, inf];
else
    edges = [-inf inf];
end
sec = zeros(nr,nc);
sec(foot) = discretize(along,edges);   % FIND and logical indexing share the same (column-major) order
sec(isnan(sec)) = 0;

DX = HAR;
DX.Z = zeros(size(HAR.Z));
DX.Z(r1:r2,c1:c2) = sec;

T = HAR;
T.Z = false(size(HAR.Z));
T.Z(r1:r2,c1:c2) = foot & sec > 0;
end
