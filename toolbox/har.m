function REM=har(DEM,S,varargin)

% HAR computing Height Above River (HAR) or Relative Elevation Model (REM)
%
% Syntax
%
%    REM=har(DEM,S)
%    REM=har(DEM,S,'power',2,'k',12,'plot',true)
%
% Description
%
%     This function allows to compute the Relative Elevation Model (REM) or
%     Height Above River (HAR) model of a river. The elevation of the river
%     course (thalweg) is sampled from the DEM and interpolated over the
%     whole grid with the Inverse Distance Weighting (IDW) of the k nearest
%     river points. The REM is the difference between the DEM and the
%     interpolated surface.
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     S         river course, as a STREAMobj or as a structure array obtained by a
%               shapefile
%
% Optional Input arguments
%
%     power     IDW power parameter (default 2)
%     k         number of nearest river points used by the IDW (default 12)
%     plot      plotting the results (default true)
%
% Output arguments
%
%     REM       high resolution GRIDobj
%
% Note: the z values of the river course are sampled from the DEM with
% INTERP, so the STREAMobj does not need to be derived from the same grid
% as the DEM. If S is a STREAMobj, z is sampled at every node (one per
% cell). If S is a shapefile, z is sampled ONLY at the vertices of the line,
% so the result depends on how the line has been digitized. Advice:
%   - draw the line over the water surface (or the channel bed), never
%     over the banks or the adjacent land;
%   - use regularly and densely spaced vertices, e.g. one vertex every
%     channel width or less (in QGIS: Processing Toolbox > Densify by
%     interval). Irregular vertices (dense in bends, sparse in straight
%     reaches) bias the IDW toward the more densely sampled reaches;
%   - use a single line feature, oriented from upstream to downstream
%     (in QGIS: Dissolve / Merge lines), especially if the line comes
%     from OpenStreetMap, where rivers are split into many features.
%
% See also: GRIDobj/interp, knnsearch
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023, - updated 27. November, 2026

% Parse inputs
p = inputParser;
p.FunctionName = 'har';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj') || isstruct(x));
addParameter(p,'power',2,@(x) isnumeric(x) && isscalar(x) && x>0)
addParameter(p,'k',12,@(x) isnumeric(x) && isscalar(x) && x>=1)
addParameter(p,'plot',true)
parse(p,DEM,S,varargin{:});

pw  = p.Results.power;
k   = p.Results.k;
plt = p.Results.plot;

%% x, y, z coordinates of the river course (thalweg)
if isa(S,'STREAMobj')
    trunk_river=trunk(klargestconncomps(S,1)); % extract trunk river
    [x,y] = STREAMobj2XY(trunk_river);          % x and y coordinates of the trunk river
elseif isstruct(S)
    x=[S.X]; % extract x coordinate
    y=[S.Y]; % extract y coordinate
end

x = double(x(:));
y = double(y(:));
z = nan(size(x));
v = ~isnan(x) & ~isnan(y);
z(v) = double(interp(DEM,x(v),y(v)));  % z sampled from the DEM

% Consistent removal of NaN values
I = ~(isnan(x) | isnan(y) | isnan(z));
x = x(I);
y = y(I);
z = z(I);

if isempty(z)
    error('har:nopoints','The river course does not overlap the DEM.')
end

k = min(round(k),numel(z));

%% IDW interpolation of the thalweg elevation
P = [x y];   % thalweg points (N x 2)

[X,Y] = getcoordinates(DEM);
[Xg,Yg] = meshgrid(X(:)', Y(:));

iq = find(~isnan(DEM.Z));      % interpolate only valid DEM cells
Zidw = nan(size(DEM.Z));

NS = KDTreeSearcher(P);
chunk = 1e6;                   % process the grid in blocks to limit memory
for c = 1:chunk:numel(iq)
    ii = iq(c:min(c+chunk-1,numel(iq)));
    [idx,D] = knnsearch(NS,[Xg(ii) Yg(ii)],'K',k);
    D(D == 0) = eps;           % avoid zero distances
    W = 1 ./ (D .^ pw);        % IDW weights
    W = W ./ sum(W,2);
    Zidw(ii) = sum(W .* reshape(z(idx),size(idx)),2);
end

% interpolated surface
A = DEM;
A.Z = cast(Zidw,'like',DEM.Z);

%% Relative Elevation Model (REM)
REM = DEM - A;

%% plotting
if plt
    figure
    imageschs(DEM, REM)
    hold on
    if isa(S,'STREAMobj')
        plot(trunk_river,'k')
    else
        plot(x,y,'k')
    end
    hc = colorbar;
    hc.Label.String = "Elevation above the thalweg (m)";
end

end
