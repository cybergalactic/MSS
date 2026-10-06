function proj = projectToPath(w_path, x_path, y_path, x, y, pi_h, k_old, dsFwd)
% projectToPath computes the orthogonal projection of a vehicle position
% onto a sampled 2-D path. The path is represented by discrete North-East
% coordinates and the corresponding path tangent angle at each sample.
%
% The search is local: only the current segment index k_old and a limited
% number of forward samples dsFwd are searched. This prevents the
% projection from jumping to a distant part of the path when the path
% crosses itself or loops.
%
% Inputs:
%   w_path   : Path parameter vector (e.g. sample index or arc-length-like
%              parameter) corresponding to x_path and y_path.
%
%   x_path   : North coordinates of the sampled path points.
%
%   y_path   : East coordinates of the sampled path points.
%
%   x, y     : Vehicle position in North-East coordinates:
%                  x = North position (m)
%                  y = East position  (m)
%
%   pi_h     : Path tangent angle at each sampled path point (rad),
%              measured relative to North:
%                  pi_h = atan2(dy_path, dx_path)
%
%   k_old    : Previous projection segment index. The search starts near
%              this segment to ensure continuity of the projection.
%
%   dsFwd    : Number of forward path samples to search when computing the
%              minimum cross-track distance.
%
% Output:
%   proj     : Struct containing the local path projection information:
%
%      proj.k
%           Segment index k such that the projection lies on the segment
%           between samples k and k+1.
%
%      proj.mu
%           Local interpolation parameter on the segment:
%               mu = 0   -> projection at sample k
%               mu = 1   -> projection at sample k+1
%
%      proj.s
%           Interpolated path parameter at the projected point.
%
%      proj.xp, proj.yp
%           North-East coordinates of the projected point on the path (m).
%
%      proj.pi
%           Path tangent angle at the projected point (rad).
%
%      proj.ex
%           Along-track error (m), positive in the path tangent direction.
%
%      proj.ey
%           Signed cross-track error (m), positive to the left of the path
%           tangent direction.
%
%      proj.dist2
%           Squared Euclidean distance from the vehicle position to the
%           projected point:
%
%               dist2 = ex^2 + ey^2
%
% Example:
%   proj = projectToPath(w_path, x_path, y_path, xN, yE, pi_h, k_old, 200);
%
% Author:    Thor I. Fossen
% Date:      2026-03-30

    backSteps = 1; % Backward tolerance in samples (default 1)

    px = x_path(:);
    py = y_path(:);
    ps = w_path(:);

    ppi = pi_h;

    N = length(ps);
    if N < 2
        error('Path must have at least two samples.');
    end

    k_old = max(1, min(k_old, N-1));
    sMax  = ps(k_old) + dsFwd;

    kMin = max(1, k_old - backSteps);
    kMax = find(ps <= sMax, 1, 'last');
    if isempty(kMax)
        kMax = min(N-1, k_old+1);
    else
        kMax = min(N-1, max(kMax, k_old+1));
    end

    p = [x; y];

    bestDist2 = inf;
    best = struct('k',k_old,'mu',0,'s',ps(k_old), ...
                  'xp',px(k_old),'yp',py(k_old), ...
                  'pi',ppi(k_old),'ex',0,'ey',0,'dist2',inf);

    for k = kMin:kMax
        pk  = [px(k);   py(k)];
        pk1 = [px(k+1); py(k+1)];
        dk  = pk1 - pk;
        L2  = dk.'*dk;

        if L2 < 1e-12
            mu = 0;
        else
            mu = ((p - pk).' * dk) / L2;
            mu = min(max(mu,0),1);
        end

        pproj = pk + mu*dk;
        d = p - pproj;
        dist2 = d.'*d;

        if dist2 < bestDist2
            sproj = ps(k) + mu*(ps(k+1)-ps(k));
            piproj = ppi(k) + mu*ssa(ppi(k+1)-ppi(k));

            t = [cos(piproj); sin(piproj)];
            n = [-sin(piproj); cos(piproj)];

            ex = t.' * d;
            ey = n.' * d;

            bestDist2 = dist2;
            best.k = k;
            best.mu = mu;
            best.s = sproj;
            best.xp = pproj(1);
            best.yp = pproj(2);
            best.pi = piproj;
            best.ex = ex;
            best.ey = ey;
            best.dist2 = dist2;
        end
    end

    proj = best;
end
