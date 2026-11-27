function bank=bankfull(d,z,varargin)

% BANKFULL extracting flow height above thalveg vs. hydraulic depth values
% from the profiles taken from the river course. 
% 
% Syntax
%
%    bank=bankfull(d,z)
%
% Description
%
%     This function allows to extract flow height above thalveg vs. 
%     hydraulic depth values from the profiles obtained by the function PROF.
%     For each profile, the hydraulic depth, defined as the ratio
%     between the flow area and the flow width, is computed for different
%     heights from the thalveg. In detail, the heigths range from 0.1 to
%     max_depth every 10 cm. Flow area and width are computed exactly on
%     the vertices of the profile and on its intersections with the water
%     level (no resampling of the profile).
%
%
% Input arguments
%
%     d         an array of the distance of the profiles
%     z         z-values of the profiles (from HAR)
%
% Optional Input arguments
%
%     max_depth maximum height from the thalveg for which the routine is
%               computed
%     connected true (default): only the wet portion of the profile
%               containing the thalweg is considered. false: all the
%               portions of the profile below the water level are
%               considered (original behaviour, it includes ditches,
%               ponds, oxbows or the other limb of a meander)
%     thalwegwindow half-width in meters of the window around the centre
%               of the profile (d = 0) where the thalweg is searched as the
%               lowest point (default 10). The wet area connected to the
%               thalweg is considered when 'connected' is true. Increase it
%               when the river course does not follow the channel, e.g.
%               when profiles of a single river course are used to compare
%               DEMs of different dates after channel migration (the window
%               should cover the channel shift, but not include gravel pits
%               or other depressions far from the channel)
%     plot      plotting the results
%
%     Each profile is closed at both ends with vertical walls higher than
%     max_depth, so the water is always confined within the profile. The
%     spill elevation of each profile (zspill, the lowest of the two highest
%     points on the left and on the right of the thalweg) is returned:
%     if zspill < max_depth the profile does not contain the valley sides
%     up to max_depth and the levels above zspill are confined by the
%     artificial walls. In this case consider increasing the profile width.
%
% Output arguments
%
%     bank      the structure array with elevation above thalveg values
%               (elevation), hydraulic depth values (h), area
%               values (area), and width values (width) for
%               each profile, spill elevation of each profile (zspill),
%               position (dthalweg) and elevation (zthalweg) of the thalweg
%
% 
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023, - updated 27. November, 2026
%
% Parse inputs
p = inputParser;
p.FunctionName = 'bankfull';
addRequired(p,'d',@(x) isnumeric(x))
addRequired(p,'z',@(x) isnumeric(x))
addParameter(p,'max_depth',6,@(x) isnumeric(x))
addParameter(p,'connected',true,@(x) isscalar(x))
addParameter(p,'thalwegwindow',10,@(x) isnumeric(x) && isscalar(x) && x>=0)
addParameter(p,'plot',true);

parse(p,d,z,varargin{:});

%% flow height above thalveg vs. hydraulic depth computation

plt=p.Results.plot;
max_depth=p.Results.max_depth; 
connected=logical(p.Results.connected);
twin=p.Results.thalwegwindow;
d_elev=0.1:0.1:max_depth;

sz=size(z);

% initialize the output variables

nlev=numel(d_elev);
h=zeros(nlev,sz(1,1));
a=zeros(nlev,sz(1,1));
w=zeros(nlev,sz(1,1));
zspill=nan(1,sz(1,1));
zthalweg=nan(1,sz(1,1));
dthalweg=nan(1,sz(1,1));

verbose = waitbar(0, 'Please wait...'); % initialize the waitbar

tic

for p=1:sz(1,1)

    y1 = z(p,:);

    % valid points of the profile (NaN: clipped ends or nodata)
    v = ~isnan(d) & ~isnan(y1);
    if nnz(v) < 2
        h(:,p)=NaN; a(:,p)=NaN; w(:,p)=NaN;
        continue
    end
    xv = d(v);
    yv = y1(v);

    % thalweg: lowest point within the window around the centre
    xt = thalwegposition(xv,yv,twin);
    dthalweg(p) = xt;
    zthalweg(p) = yv(find(xv==xt,1));

    % spill elevation of the profile (diagnostic)
    zl = max(yv(xv<=xt)); zr = max(yv(xv>=xt));
    if ~isempty(zl) && ~isempty(zr)
        zspill(p) = min(zl,zr);
    end

    % flow area and width for each level (exact, profile closed by walls)
    for i=1:nlev
        [area,width] = wetgeometry(xv,yv,d_elev(i),connected,xt);
        h(i,p)=area./width;
        a(i,p)=area;
        w(i,p)=width;
    end

    waitbar(p/sz(1,1), verbose, sprintf('Processing section %d of %d', p, sz(1,1)));
end

close(verbose)

toc



bank.elevation=d_elev;
bank.h=h;
bank.area=a;
bank.width=w;
bank.zspill=zspill;
bank.dthalweg=dthalweg;
bank.zthalweg=zthalweg;

bank.area(isinf(bank.h))=nan;
bank.width(isinf(bank.h))=nan;
bank.h(isinf(bank.h))=nan;

% warn if the thalweg found is not on the river
noff = sum(zthalweg > 1);
if noff>0
    warning('bankfull:thalweg',['In %d of %d profiles the thalweg (lowest point within %g m from the centre) is ' ...
        'more than 1 m above the river course (median %.2f m): the profiles may not be centred on the channel. ' ...
        'Check the river course, use ''keeptrace'',true in PROF if the trace is smoothed, or increase ''thalwegwindow''.'], ...
        noff,sz(1,1),twin,median(zthalweg(zthalweg>1)));
end

% warn if the profiles do not contain the valley sides up to max_depth
nshort=sum(zspill<max_depth);
if nshort>0
    warning('bankfull:spill',['%d of %d profiles have a spill elevation (zspill) lower than max_depth (%.1f m): ' ...
        'above zspill the water is confined only by the artificial walls closing the profile. ' ...
        'Consider increasing the profile width in PROF (median zspill = %.2f m).'], ...
        nshort,sz(1,1),max_depth,median(zspill,'omitnan'));
end

%% plotting results

if plt

c=parula(sz(1,1));

figure

for p=1:sz(1,1)

    plot(d_elev,h(:,p),'Color',c(p,:))
    hold on
    scatter(d_elev,h(:,p),20,'MarkerFaceColor',c(p,:),'MarkerEdgeColor','none','Marker','o','SizeData',10,'LineWidth',.1)

end

set(gca, 'YScale', 'log')
ylabel('Hydraulic depth (area/width) (m)')
xlabel ('Elevation above thalweg (m)')

end

end

%% ------------------------------------------------------------------------
function [A,W,P] = wetgeometry(xv,yv,L,connected,xt)
% WETGEOMETRY exact flow area (A), flow width (W) and wetted perimeter (P)
% of the profile (xv,yv) for the water level L.
%
% The profile is a polyline (the bed is linear between consecutive
% vertices) closed at both ends by vertical walls, so the water is always
% confined within the profile. Each wet interval is delimited by the
% intersections between the water level and the profile (or by the walls);
% the area is computed with TRAPZ on the vertices and the intersections,
% which is exact for a polyline. The walls are not part of the wetted
% perimeter. If connected is true, only the wet interval containing the
% thalweg (x = xt) is considered.

A = 0; W = 0; P = 0;

[xv,o] = sort(xv(:)');
yv = yv(:)';
yv = yv(o);
n = numel(xv);

wet = yv < L;
if ~any(wet)
    return
end

% runs of consecutive wet vertices
dw = diff([false wet false]);
i1 = find(dw == 1);
i2 = find(dw == -1) - 1;

for k = 1:numel(i1)
    a = i1(k);
    b = i2(k);
    xs = xv(a:b);
    ys = yv(a:b);

    % left boundary: intersection with the water level or wall
    if a > 1
        xa = xv(a-1) + (L-yv(a-1)) / (yv(a)-yv(a-1)) * (xv(a)-xv(a-1));
        xs = [xa xs]; %#ok<AGROW>
        ys = [L ys];  %#ok<AGROW>
    end
    % right boundary: intersection with the water level or wall
    if b < n
        xb = xv(b) + (L-yv(b)) / (yv(b+1)-yv(b)) * (xv(b+1)-xv(b));
        xs = [xs xb]; %#ok<AGROW>
        ys = [ys L];  %#ok<AGROW>
    end

    if connected && ~(xs(1) <= xt && xs(end) >= xt)
        continue
    end

    A = A + trapz(xs, L - ys);
    W = W + (xs(end) - xs(1));
    P = P + sum(hypot(diff(xs), diff(ys)));

    if connected
        break
    end
end
end

%% ------------------------------------------------------------------------
function xt = thalwegposition(xv,yv,twin)
% THALWEGPOSITION position of the thalweg along the profile (xv,yv): the
% lowest point within twin meters from the centre of the profile (d = 0).
% If no point lies within the window, the point closest to the centre is
% used.
m = abs(xv) <= twin;
if ~any(m)
    [~,k] = min(abs(xv));
    xt = xv(k);
    return
end
xm = xv(m);
ym = yv(m);
[~,k] = min(ym);
xt = xm(k);
end
