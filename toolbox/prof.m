function [d,z,z_dem,x,y,SW]=prof(DEM,HAR,S,varargin)

% PROF extracting profiles from SWATHobj of the river course
%
% Syntax
%
%    [d,z,z_dem,x,y,SW]=prof(DEM,HAR,S)
%
% Description
%
%     This function allows to extract profile from Relative Elevation
%     Model (REM) or Height Above River (HAR) model as well as from Digital
%     Elevation Model (DEM) 
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     HAR       high resolution GRIDobj
%     S         river course, as a STREAMobj or as a structure array obtained by a
%               shapefile or from an free hand drawn line
%
% Optional Input arguments
%     step      stepping in meters between transversal river profiles
%     width     width in meters of the profiles
%     smooth    smoothing of the trace from which the profiles are
%               extracted: false (default) or true. In TopoToolbox 3 the
%               trace is smoothed with SMOOTHDATA (see SWATHobj); the
%               strength is set by 'smoothingfactor'. A numeric value is
%               interpreted as true if > 0 (old syntax, in meters, is no
%               longer supported by TopoToolbox 3).
%     smoothingfactor  value between 0 and 1 controlling the smoothing of
%               the trace (default 0.25, see SMOOTHDATA)
%     keeptrace true (default): if the trace is smoothed, the smoothed
%               trace is used only to orient the transverse profiles, which
%               remain centred on the original river course (d = 0 on the
%               river). false: the profiles are centred on the smoothed
%               trace, which can depart from the channel in sinuous reaches
%     clip      true (default): each transverse profile is cut where it
%               enters the "territory" of another reach of the river, i.e.
%               where its points are closer to a point of the river course
%               located more than 'cliptol' meters away (along the course)
%               from the profile centre. This prevents profiles from
%               crossing the other limb of tight meanders. The removed
%               points are set to NaN. false: original behaviour
%     cliptol   along-course tolerance in meters for 'clip' (default: the
%               profile width)
%     plot      plotting the results
%
% Output arguments
%
%     d         an array of the distance of the profiles
%     z         z-values of the profiles (from HAR)
%     z_dem     z-values of the profiles (from DEM)
%     x         x-coordinates of the profiles
%     y         y-coordinates of the profiles
%     SW        SWATHobj of the river profile
%
% See also: SWATHobj/profiles
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023, - updated 27. November, 2026

% Parse inputs
p = inputParser;
p.FunctionName = 'prof';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'HAR',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj') || isstruct(x)||ischar(validatestring(x,{'interactive'})));
addParameter(p,'step',1e3,@(x) isnumeric(x))
addParameter(p,'width',1e3,@(x) isnumeric(x))
addParameter(p,'smooth',false,@(x) isscalar(x) && (islogical(x) || isnumeric(x)))
addParameter(p,'smoothingfactor',0.25,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1)
addParameter(p,'keeptrace',true,@(x) isscalar(x))
addParameter(p,'clip',true,@(x) isscalar(x))
addParameter(p,'cliptol',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0))
addParameter(p,'plot',true);

parse(p,DEM,HAR,S,varargin{:});


S=p.Results.S;
plt=p.Results.plot;
step=p.Results.step;

% smoothing of the trace (TopoToolbox 3: logical switch + smoothing factor)
smooth=p.Results.smooth;
if isnumeric(smooth) && smooth>1
    warning('prof:smooth',['In TopoToolbox 3 ''smooth'' is a switch (true/false), not a distance. ' ...
        'smooth=%g is interpreted as true; use ''smoothingfactor'' (0-1) to control the smoothing.'],smooth)
end
smooth=logical(smooth);
sf=p.Results.smoothingfactor;
keeptrace=logical(p.Results.keeptrace);

clip=logical(p.Results.clip);
cliptol=p.Results.cliptol;
if isempty(cliptol)
    cliptol=p.Results.width;
end
%% extracting profiles

if isa(S,'STREAMobj')

    [X,Y] = STREAMobj2XY(S); % extract x and y coordinates of the river
    
    ix = ~isnan(X); % create a index array to remove nans from the coordinates

    SW = SWATHobj(DEM,X(ix),Y(ix),'width',p.Results.width,'smooth',smooth,'smoothingfactor',sf,'keeptrace',keeptrace); % compute the swath along the river

    if clip
        SW = cliptransversals(SW,cliptol); % cut the profiles at other reaches of the river
    end

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath

    SW = mapswath(SW,HAR); % map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath

elseif isstruct(S)

    X=[S.X]; % extract x coordinate of the river

    Y=[S.Y]; % extract y coordinate of the river

    ix = ~isnan(X); % create a index array to remove nans from the coordinates

    SW = SWATHobj(DEM,X(ix),Y(ix),'width',p.Results.width,'smooth',smooth,'smoothingfactor',sf,'keeptrace',keeptrace); % compute the swath along the river

    if clip
        SW = cliptransversals(SW,cliptol); % cut the profiles at other reaches of the river
    end

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath

    SW = mapswath(SW,HAR); % map swath with HAR% map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath
    

elseif ischar(S)

    SW= SWATHobj(DEM,'width',p.Results.width); % interactively digitize the river and compute the swath along the river

    if clip
        SW = cliptransversals(SW,cliptol); % cut the profiles at other reaches of the river
    end

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath  

    SW = mapswath(SW,HAR); % map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath

end

%% plotting profiles
if plt

figure

imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false)
hold on

c=parula(size(z_dem,1));

X = reshape(x,[],size(z_dem,1));
Y = reshape(y,[],size(z_dem,1));

for i =1:size(z_dem,1)
    plot(X(:,i),Y(:,i),'Color',c(i,:))
    hold on
end

plot(SW)

% xlim([3.467e5 3.615e5])
% ylim([4.786e6 4.7925e6])

figure

for i =1:size(z_dem,1)
    plot(d,z_dem(i,:),'Color',c(i,:))
    hold on
end

ylabel('Elevation (m a.s.l.)')
xlabel('Distance from thalweg (m)')

end

end

%% ------------------------------------------------------------------------
function SW = cliptransversals(SW,tol)
% Cut each transverse profile of the SWATHobj where its points are closer
% to a point of the trace located more than tol meters away (along the
% trace) from the profile centre. Moving outward from the centre, all the
% points beyond the first such point are set to NaN, so that each profile
% remains a continuous segment centred on the trace.

[ny,nx] = size(SW.Z);
dist = SW.distx(:);

% nearest trace point of each point of the swath
idx = knnsearch(SW.xy,[SW.X(:) SW.Y(:)]);
own = repmat(1:nx,ny,1);
bad = abs(dist(idx) - dist(own(:))) > tol;
bad = reshape(bad,ny,nx);

% centre of the transverse profiles
[~,ic] = min(abs(SW.disty));

cut = false(ny,nx);
for j = 1:nx
    % left side (from the centre outward)
    k = find(bad(ic:-1:1,j),1,'first');
    if ~isempty(k)
        cut(1:ic-k+1,j) = true;
    end
    % right side (from the centre outward)
    k = find(bad(ic:ny,j),1,'first');
    if ~isempty(k)
        cut(ic+k-1:ny,j) = true;
    end
end

SW.Z(cut) = NaN;
end
