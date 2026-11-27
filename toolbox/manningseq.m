function R=manningseq(SW,z,z_dem,step,d,zvar,varargin)

% MANNINGSEQ computes the Manning's Equation Q = 1/n * (A/P)^(2/3) * S^(1/2)
% * A
%
% Syntax
%
%    R=manningseq(SW,z,z_dem,step,d,H)
%
% Description
%
%     This function allows to compute the river discharge using the
%     Manning's equation. Flow area, width and wetted perimeter are
%     computed exactly on the vertices of the profile and on its
%     intersections with the water level (no resampling of the profile).
%
% Input arguments
%
%     SW        SWATHobj of the river profile
%     z         z-values of the profiles (from HAR)
%     z_dem     z-values of the profiles (from DEM)
%     step      distance in meters between the profiles
%     d         an array of the distance of the profiles
%     zvar      variable height above thalveg (from DETECT_PEAK function)
%     slope     river gradient
%
% Optional Input arguments
%     n         n coefficient
%     slope     river gradient: scalar or one value per section
%     connected true (default): flow area, width and wetted perimeter are
%               computed only on the wet portion containing the thalweg.
%               false: original behaviour (all portions below zvar)
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
% Output arguments
%
% The outputs are organized in a struct file where are organized as follows:
%
%     Q_zvar     river discharge for each section
%     Qi_zvar    most probable river discharge
%     area_zvar  flow area for each section
%     wet_p_zvar wet perimeter for each section
%     width_zvar width for each section
%     slope      river gradient for each section
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024, - updated 27. November, 2026
%
%
% Parse inputs

p = inputParser;
p.FunctionName = 'manningseq';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addRequired(p,'z',@(x) isnumeric(x))
addRequired(p,'z_dem',@(x) isnumeric(x))
addRequired(p,'step',@(x) isnumeric(x))
addRequired(p,'d',@(x) isnumeric(x))
addRequired(p,'zvar',@(x) isnumeric(x))
addParameter(p,'n',0.05,@(x) isnumeric(x))
addParameter(p,'slope',[],@(x) isnumeric(x));
addParameter(p,'connected',true,@(x) isscalar(x))
addParameter(p,'thalwegwindow',10,@(x) isnumeric(x) && isscalar(x) && x>=0)
addParameter(p,'plot',true)

parse(p,SW,z,z_dem,step,d,zvar,varargin{:});

n=p.Results.n;
slope=p.Results.slope;
plt=p.Results.plot;
connected=logical(p.Results.connected);
twin=p.Results.thalwegwindow;
zvar=zvar(:)';   % row vector

%%

sz=size(z);

% initialize the variables

area_zvar=zeros(1,sz(1,1));
wet_p_zvar=zeros(1,sz(1,1));
width_zvar=zeros(1,sz(1,1));

verbose = waitbar(0, 'Please wait...'); % initialize the waitbar

tic

for p=1:sz(1,1)

    y1 = z(p,:);

    % valid points of the profile (NaN: clipped ends or nodata)
    v = ~isnan(d) & ~isnan(y1);

    if isnan(zvar(1,p)) || nnz(v) < 2
        % Skip when zvar(1,p) is NaN or the profile is empty
        area_zvar(1,p)=NaN;
        wet_p_zvar(1,p)=NaN;
        width_zvar(1,p)=NaN;
        continue
    end

    % flow area, width and wetted perimeter (exact, profile closed by walls)
    xt = thalwegposition(d(v),y1(v),twin);
    [area_zvar(1,p),width_zvar(1,p),wet_p_zvar(1,p)] = wetgeometry(d(v),y1(v),zvar(1,p),connected,xt);

    waitbar(p/sz(1,1), verbose, sprintf('Processing section %d of %d', p, sz(1,1)));

end

close(verbose)

toc

if isempty(slope)

    slope=abs((SW.zd0(end,1)-SW.zd0(1,1))/SW.zd0(end,2));

end

slope=slope(:)';   % row vector (one value per section) or scalar
if ~isscalar(slope) && numel(slope)~=sz(1,1)
    error('manningseq:slope','slope must be a scalar or have one value per section (%d values given, %d sections).',numel(slope),sz(1,1))
end

Q_zvar=1/n.*((area_zvar./wet_p_zvar).^(2/3)).*(abs(slope).^(1/2)).*area_zvar;

Q_zvar(Q_zvar==0)=NaN;

[f, xi] = ksdensity(Q_zvar(~isnan(Q_zvar)), 'Bandwidth', 10);

Qi_zvar=xi(find(f==max(f),1));

% make struct file of the results

R.Q_zvar=Q_zvar;
R.Qi_zvar=Qi_zvar;

R.area_zvar=area_zvar;
R.wet_p_zvar=wet_p_zvar;

R.width_zvar=width_zvar;

R.slope=slope;

if plt

    % plotting histogram Q
    figure

    histogram(Q_zvar,100)

    ylabel('Count')
    xlabel('River Discharge (m^{3} s^{-1})')

    yyaxis right

    plot(xi,f,'r-','LineWidth',1)
    ylabel ('Density probability')


    xline(Qi_zvar,':','LineWidth',1);
    

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
