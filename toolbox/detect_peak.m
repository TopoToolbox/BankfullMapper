function lim=detect_peak(h,a,w,d_elev,d,z,varargin)

% DETECT_PEAK detecting peaks from the flow height above thalveg vs.
% hydraulic depth curves of each profile
%
% Syntax
%
%    lim=detect_peak(h,a,w,d_elev,d,z)
%
% Description
%
%     This function allows to detect the peaks from the flow heigh above thalveg vs. hydraulic depth
%     hydraulic depth curves of each profiles. It is possible to detect the
%     all the peaks (using 'all' option), most prominent peaks (using 'max' option),
%     and lowest peaks (using 'lowest' option). For the 'max' and 'lowest'
%     option, the most probable Z value is also computed.
%
%
% Input arguments
%
%     h         hydraulic depth values
%     a         flow area values
%     w         flow width values
%     d_elev    flow height values
%
% Optional Input arguments
%     d         an array of the distance of the profiles (from prof
%               function). It is necessary if 'peak' option is 'max' or 'lowest'
%     z         z-values of the profiles (from HAR). It is necessary if 'peak' option is 'max' or 'lowest'
%     'peak'    string
%        'all': all the peaks are extracted
%        'max': the most prominent peaks are extracted
%        'lowest': the lowest peaks are extracted
%     connected  true (default): flow area and width at zvar are computed
%                only on the wet portion containing the thalweg.
%                false: original behaviour (all portions below zvar)
%     thalwegwindow half-width in meters of the window around the centre
%               of the profile (d = 0) where the thalweg is searched as the
%               lowest point (default 10). The wet area connected to the
%               thalweg is considered when 'connected' is true. Increase it
%               when the river course does not follow the channel, e.g.
%               when profiles of a single river course are used to compare
%               DEMs of different dates after channel migration (the window
%               should cover the channel shift, but not include gravel pits
%               or other depressions far from the channel)
%      plot      plotting the results
%
%     Peaks are detected with ISLOCALMAX (MATLAB base): the Signal
%     Processing Toolbox is not required. Only the interior local maxima of
%     the first continuous run of non-NaN values of each curve are
%     considered. With the 'max' option, if two peaks have the same
%     prominence the lowest one is taken.
%
% Output arguments
%
%     lim       the structure array with section numbers (sections field),
%               peak elevation above thalveg values (zlim), the most probable 
%               peak elevation above thalveg for each section based on the 
%               2d kernel distribution estimation of zlim values (zvar), peak hydraulic
%               depth values (hlim), peak hydraulic depth values using zvar 
%               (hvar only if 'peak' option is 'max' or 'lowest'), flow area values (alim),
%               flow area values using zvar (avar only if 'peak' option is 'max' or 'lowest'), 
%               flow width values (wlim), flow width values using zvar (wvar only if 'peak' option is 'max' or 'lowest'). 
%               For 'max' and 'lowest' peak detection, the most probable elevation above thalveg
%               value of the peaks is also reported (Z).
%
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024, - updated 27. November, 2026
%
% Parse inputs
p = inputParser;
p.FunctionName = 'detect_peak';

addRequired(p,'h',@(x) isnumeric(x))
addRequired(p,'a',@(x) isnumeric(x))
addRequired(p,'w',@(x) isnumeric(x))
addRequired(p,'d_elev',@(x) isnumeric(x))
addOptional(p,'d',@(x) isnumeric(x))
addOptional(p,'z',@(x) isnumeric(x))
addParameter(p,'peak','all',@(x) ischar(validatestring(x,{'all','max','lowest'})))
addParameter(p,'connected',true,@(x) isscalar(x))
addParameter(p,'thalwegwindow',10,@(x) isnumeric(x) && isscalar(x) && x>=0)
addParameter(p,'plot',true)

parse(p,h,a,w,d_elev,d,z,varargin{:});

%% peak detection
plt=p.Results.plot;
connected=logical(p.Results.connected);
twin=p.Results.thalwegwindow;

meth = validatestring(p.Results.peak,{'all','max','lowest'});

sz=size(h);

switch meth

    case 'all'

        % detect zlim and hlim values

        zl = cell(sz(1,2), 1);
        hl = cell(sz(1,2), 1);


        % Store the peak elevation and hydraulic depth values in zlim
        % and hlim vectors

        zlim = [];
        hlim = [];
        alim = [];
        wlim = [];

        % Store the sections numbers in the sections vector
        sections = [];

        for p = 1:sz(1,2)

            idx = localpeaks(h(:,p));
            hl{p,1} = h(idx,p);
            zl{p,1} = d_elev(idx);

            % Concatenate the values to the zlim and hlim arrays
            zlim = [zlim, zl{p, 1}];
            hlim = [hlim; hl{p, 1}];
            alim = [alim, a(idx,p)'];
            wlim = [wlim, w(idx,p)'];

            % Store the sections numbers
            sections = [sections, repmat(p, 1, numel(zl{p, 1}))];


        end

        % remove infinite values (if any)
        ii = isinf(hlim(:)');
        zlim(ii)=[];
        sections(ii)=[];
        alim(ii)=[];
        wlim(ii)=[];
        hlim(ii)=[];


        lim.sections=sections;
        lim.zlim=zlim;
        lim.hlim=hlim';
        lim.alim=alim;
        lim.wlim=wlim;



    case 'max'

        sections=1:1:sz(1,2);
        zlim=zeros(1,sz(1,2));
        hlim=zeros(1,sz(1,2));
        h_var=zeros(1,sz(1,2));
        alim=zeros(1,sz(1,2));
        area_zvar=zeros(1,sz(1,2));
        width_zvar=zeros(1,sz(1,2));
        wlim=zeros(1,sz(1,2));



        for p=1:sz(1,2)
            wt=w(:,p);
            at=a(:,p);

            % Find the indices and prominence of the local maxima
            [idx,pr] = localpeaks(h(:,p));

            if isempty(idx)
                zlim(1,p)=nan;
                hlim(1,p)=nan;
                alim(1,p)=NaN;
                wlim(1,p)=NaN;
            else
                % most prominent peak (the lowest one if more peaks have
                % the same prominence)
                [~,j] = max(pr);
                k = idx(j);
                zlim(1,p)=d_elev(k);
                hlim(1,p)=h(k,p);
                alim(1,p)=at(k);
                wlim(1,p)=wt(k);
            end


        end

        [f, xi] = ksdensity(zlim(~isnan(zlim)), 'Bandwidth', 0.2);
        Z=xi(find(f==max(f),1));

        % computing the zvar from the most probable kernel distribution
        n_points = 100;
        [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
        xy = [sections', zlim'];
        xy = xy(~any(isnan(xy),2),:);   % sections without peak are ignored
        [f,~] = ksdensity(xy, [xi(:), yi(:)]);
        f = reshape(f, n_points, n_points);

        % Initialize array to store mode zvar for each section
        zvar = zeros(n_points, 1);

        for i = 1:n_points
            % Find the index of the maximum density value for the current section
            [~, max_idx] = max(f(:, i));

            % Get the corresponding zvar value
            zvar(i) = yi(max_idx, i);
        end

        % Interpolate mode_zvar values for the full length of sections
        zvar = interp1(linspace(min(sections), max(sections), n_points), zvar, sections, 'linear');

        %%
        for p=1:sz(1,2)

            y1 = z(p,:);

            % valid points of the profile (NaN: clipped ends or nodata)
            v = ~isnan(d) & ~isnan(y1);

            if isnan(zvar(1,p)) || nnz(v) < 2
                width_zvar(1,p)=NaN;
                area_zvar(1,p)=NaN;
                h_var(1,p)=NaN;
                continue
            end

            % flow area and width at zvar (exact, profile closed by walls)
            xt = thalwegposition(d(v),y1(v),twin);
            [area_zvar(1,p),width_zvar(1,p)] = wetgeometry(d(v),y1(v),zvar(1,p),connected,xt);
            h_var(1,p)=area_zvar(1,p)./width_zvar(1,p);

        end

        %%
        lim.sections=sections;
        lim.zlim=zlim;
        lim.zvar=zvar;
        lim.hlim=hlim;
        lim.hvar=h_var;
        lim.alim=alim;
        lim.avar=area_zvar;
        lim.wlim=wlim;
        lim.wvar=width_zvar;
        lim.Z=Z;





    case 'lowest'

        sections=1:1:sz(1,2);
        zlim=zeros(1,sz(1,2));
        hlim=zeros(1,sz(1,2));
        h_var=zeros(1,sz(1,2));
        alim=zeros(1,sz(1,2));
        area_zvar=zeros(1,sz(1,2));
        width_zvar=zeros(1,sz(1,2));
        wlim=zeros(1,sz(1,2));

        for p=1:sz(1,2)
            wt=w(:,p);
            at=a(:,p);
            % Find the indices of the local maxima
            idx = localpeaks(h(:,p));

            if isempty(idx)
                zlim(1,p)=nan;
                hlim(1,p)=nan;
                alim(1,p)=NaN;
                wlim(1,p)=NaN;
            else
                % lowest peak
                k = idx(1);
                zlim(1,p)=d_elev(k);
                hlim(1,p)=h(k,p);
                alim(1,p)=at(k);
                wlim(1,p)=wt(k);
            end


        end

        [f, xi] = ksdensity(zlim(~isnan(zlim)), 'Bandwidth', 0.2);
        Z=xi(find(f==max(f),1));

        % computing the zvar from the most probable kernel distribution
        n_points = 400;
        [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
        xy = [sections', zlim'];
        xy = xy(~any(isnan(xy),2),:);   % sections without peak are ignored
        [f,~] = ksdensity(xy, [xi(:), yi(:)]);
        f = reshape(f, n_points, n_points);

        % Initialize array to store mode zvar for each section
        zvar = zeros(n_points, 1);

        for i = 1:n_points
            % Find the index of the maximum density value for the current section
            [~, max_idx] = max(f(:, i));

            % Get the corresponding zvar value
            zvar(i) = yi(max_idx, i);
        end

        % Interpolate mode_zvar values for the full length of sections
        zvar = interp1(linspace(min(sections), max(sections), n_points), zvar, sections, 'linear');


        %%
        for p=1:sz(1,2)

            y1 = z(p,:);

            % valid points of the profile (NaN: clipped ends or nodata)
            v = ~isnan(d) & ~isnan(y1);

            if isnan(zvar(1,p)) || nnz(v) < 2
                width_zvar(1,p)=NaN;
                area_zvar(1,p)=NaN;
                h_var(1,p)=NaN;
                continue
            end

            % flow area and width at zvar (exact, profile closed by walls)
            xt = thalwegposition(d(v),y1(v),twin);
            [area_zvar(1,p),width_zvar(1,p)] = wetgeometry(d(v),y1(v),zvar(1,p),connected,xt);
            h_var(1,p)=area_zvar(1,p)./width_zvar(1,p);

        end

        %%
        lim.sections=sections;
        lim.zlim=zlim;
        lim.zvar=zvar;
        lim.hlim=hlim;
        lim.hvar=h_var;
        lim.alim=alim;
        lim.avar=area_zvar;
        lim.wlim=wlim;
        lim.wvar=width_zvar;
        lim.Z=Z;

end

%% plotting

if plt
    switch meth

        case 'all'
            figure
            scatterhist(sections, zlim,'NBins',[sz(1, 2),round(max(zlim)/0.1)],'Location', 'SouthEast','Marker','.','Color','k','Kernel','overlay','Bandwidth',[2;.2]);
            xlabel('Section number');
            ylabel('Peak Elevation (m)');
            xlim([1 sz(1,2)])
            ylim([0 max(zlim)])

            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            xy = xy(~any(isnan(xy),2),:);   % sections without peak are ignored
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')



            subplot(4,1,2)

            scatter(sections,wlim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')


            subplot(4,1,4)

            scatter(sections,hlim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')



        case 'max'

            figure

            histogram(zlim,round(max(zlim)/0.1))

            ylabel('Count')
            xlabel('Elevation above thalveg (m)')

            yyaxis right

            [f, xi] = ksdensity(zlim(~isnan(zlim)), 'Bandwidth', 0.2);

            plot(xi,f,'r-','LineWidth',1)

            xline(Z,'k-','LineWidth',1);
            xline(median(zvar),'r-','LineWidth',1)
            xline(median(zvar)-std(zvar),'r--','LineWidth',1)
            xline(median(zvar)+std(zvar),'r--','LineWidth',1)

            ylabel ('Density probability')

            
            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            xy = xy(~any(isnan(xy),2),:);   % sections without peak are ignored
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')
            plot(sections,zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')


            subplot(4,1,2)

            scatter(sections,wlim,'k.')
            hold on
            plot(sections,width_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')
            hold on
            plot(sections,area_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')

            subplot(4,1,4)

            scatter(sections,hlim,'k.')
            hold on
            plot(sections,h_var, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')


        case 'lowest'

            figure

            histogram(zlim,round(max(zlim)/0.1))

            ylabel('Count')
            xlabel('Elevation above thalveg (m)')

            yyaxis right

            [f, xi] = ksdensity(zlim(~isnan(zlim)), 'Bandwidth', 0.2);

            plot(xi,f,'r-','LineWidth',1)

            xline(Z,'k-','LineWidth',1);
            xline(median(zvar),'r-','LineWidth',1)
            xline(median(zvar)-std(zvar),'r--','LineWidth',1)
            xline(median(zvar)+std(zvar),'r--','LineWidth',1)

            ylabel ('Density probability')

            
            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            xy = xy(~any(isnan(xy),2),:);   % sections without peak are ignored
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')
            plot(sections,zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')


            subplot(4,1,2)

            scatter(sections,wlim,'k.')
            hold on
            plot(sections,width_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')
            hold on
            plot(sections,area_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')

            subplot(4,1,4)

            scatter(sections,hlim,'k.')
            hold on
            plot(sections,h_var, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')



    end
end

end

%% ------------------------------------------------------------------------
function [idx,pr] = localpeaks(hc)
% interior local maxima (indices and prominence, ISLOCALMAX) of the first
% continuous run of non-NaN values of the hydraulic depth curve hc
hc = hc(:);
idx = zeros(0,1);
pr  = zeros(0,1);
i1 = find(~isnan(hc),1,'first');
if isempty(i1)
    return
end
i2 = find(isnan(hc(i1:end)),1,'first');
if isempty(i2)
    i2 = numel(hc);
else
    i2 = i1 + i2 - 2;
end
seg = hc(i1:i2);
if numel(seg) < 3
    return
end
[TF,P] = islocalmax(seg);
idx = find(TF) + i1 - 1;
pr  = P(TF);
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
