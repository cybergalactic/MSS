function testOrder3Order5TwoWaypoints()
% testOrder3Order5TwoWaypoints() checks the two-waypoint branch of order3.m
% and order5.m.
%
% For two waypoints (x1,y1), (x2,y2) both functions insert a third waypoint
% "in the middle" and return the polynomial table of the two-segment path.
% The inserted waypoint must be the midpoint ((x1+x2)/2, (y1+y2)/2): the
% first segment then ends there, and the path is the straight line between
% the waypoints. The file used (x2-x1)/2, which is the midpoint only when the
% first waypoint is the origin (the 3-or-more-waypoint branch uses absolute
% coordinates, and so does the first waypoint of this one).
%
% Checks, for first waypoints at the origin and away from it:
%   (1) the end point of the first segment equals the midpoint;
%   (2) the last point of the second segment equals the second waypoint;
%   (3) every sampled point of the path lies on the straight line between the
%       waypoints (distance to the chord), and the path length equals the
%       chord length.
% order3 is called with zero end velocities (mode 1, xstart = xfinal =
% ystart = yfinal = 0), order5 as in its header example.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Author:    Enio Krizman
% Date:      2 November 2026

cases = {[0 200],[0 100];          % first waypoint at the origin
         [100 300],[50 250];       % away from the origin
         [-40 120],[30 -90];
         [500 100],[-20 60]};

nfail = 0;
for fn = {'order3','order5'}
    for k = 1:size(cases,1)
        x = cases{k,1}; y = cases{k,2};
        if strcmp(fn{1},'order3')
            [tab,L] = order3(x,y,0,1,0,0,0,0);
            nx = 4; xi = 1:4; yi = 5:8; th2 = [1 2]; th1 = [0 1];
        else
            [tab,L] = order5(0,0,x,y);
            nx = 6; xi = 1:6; yi = 7:12; th2 = [0 1]; th1 = [0 1];
        end
        xm = mean(x); ym = mean(y);
        % (1) end of segment 1, (2) end of segment 2
        e1 = norm([polyval(tab(xi,1),th1(2)) - xm, polyval(tab(yi,1),th1(2)) - ym]);
        e2 = norm([polyval(tab(xi,2),th2(2)) - x(2), polyval(tab(yi,2),th2(2)) - y(2)]);
        % (3) distance to the chord and path length
        s = linspace(0,1,401);
        px = [polyval(tab(xi,1),th1(1) + s*diff(th1)), polyval(tab(xi,2),th2(1) + s*diff(th2))];
        py = [polyval(tab(yi,1),th1(1) + s*diff(th1)), polyval(tab(yi,2),th2(1) + s*diff(th2))];
        c = [x(2)-x(1); y(2)-y(1)]; c = c/norm(c);
        dist = max(abs((px - x(1))*c(2) - (py - y(1))*c(1)));
        eL = abs(L - hypot(x(2)-x(1), y(2)-y(1)));
        ok = e1 < 1e-9 && e2 < 1e-9 && dist < 1e-9 && eL < 1e-3;
        fprintf('%s x=%s y=%s: mid err %.3g, end err %.3g, off-line %.3g, length err %.3g  %s\n', ...
            fn{1}, mat2str(x), mat2str(y), e1, e2, dist, eL, passfail(ok));
        nfail = nfail + ~ok;
    end
end
if nfail > 0
    error('testOrder3Order5TwoWaypoints:failed', '%d of %d cases failed', nfail, 2*size(cases,1));
end
disp('testOrder3Order5TwoWaypoints PASSED');
end

function s = passfail(ok)
if ok, s = 'ok'; else, s = 'FAIL'; end
end
