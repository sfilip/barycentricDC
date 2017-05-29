function [rh,xk,err] = barycentricDC(fx, x, m, n, varargin)
%BARYCENTRICDC   Computes a best type (m,n) approximation to the
%       data fx on the set of points x.
%       R = BARYCENTRICDC(FX,X,M,N) returns a function handle R
%       for a rational approximation R expressed in barycentric
%       form (i.e., quotient N/D of partial fraction representations
%       that have the same poles).
%       [R, XK, ERR] = BARYCENTRICDC(FX,X,M,N) returns a set of points
%       where the error function FX-R equioscillates on X, while ERR
%       represents an estimation of the max value of FX-R over X.
%       [R,XK,ERR] = BARYCENTRICDC(FX,X,M,N,OPTION). The OPTION
%       parameter specifies if any execution information should
%       shared with the user during execution. The options are:
%           - 'log' - shows the approximation errors of each
%                     approximation of intermediary type
%                     (m_k,n_k), m_k <= m, n_k <= n computed
%                     during execution & the execution information
%                     when solving the LP problems;
%           - 'display' - plots the approximation error for each
%                         approximation of intermediary type (m_k,n_k)
%                         computed during execution.
[logFlag, displayFlag] = parseInputs(varargin{:});
startM = m - min(m,n);
startN = n - min(m,n);

tk = linspace(x(1), x(end), max(startM,startN)+1);
for ii = 1:length(tk)
    if(min(abs(x-tk(ii))) < 1e-8)
        tk(ii) = tk(ii) + 1e-5;
    end
end
err = 100;

while ( startM <= m )
    [rh,xk,err] = dcKernel(fx, x, startM, startN, tk, err, ...
                            logFlag, displayFlag);
    startM = startM + 1; startN = startN + 1;
    tk = refGen(x,xk,max(startM,startN)+1,0);
    for ii = 1:length(tk)
        if(min(abs(x-tk(ii))) < 1e-8)
            tk(ii) = tk(ii) + 1e-5;
        end
    end
end

end

function [logFlag, displayFlag] = parseInputs(varargin)

    logFlag = false;
    displayFlag = false;
    for k = 1:length(varargin)
        if ( strcmpi('display', varargin{k}) )
        displayFlag = true;
        elseif ( strcmpi('log', varargin{k}) )
        logFlag = true;
        else
            error('BARYDC:badInput', ...
            'Unrecognized sequence of input parameters.')
        end
    end
end



function [rh, xk, err] = dcKernel(fx, x, m, n, tk, errPrev, ...
                                  logFlag, displayFlag)

if m ~= n                               % force coefficients
                                        % to lie in null space
                                        % of Vandermonde
    mndiff = abs(m-n);
    Q = ones(length(tk),1);
    Q = Q/norm(Q);
    for ii = 2:mndiff
        Qtmp = diag(tk)*Q(:,end);
        Qtmp = Qtmp-Q*(Q'*Qtmp);        % orthogonalize
        Qtmp = Qtmp/norm(Qtmp);         % normalize
        Q = [Q Qtmp];
    end
    [Q,~] = qr(Q);
    Qmn = Q(:,1:mndiff);
end

TpowersDenom = [];

node = @(z) prod(z-tk);                 % needed to check sign
nodevec = zeros(length(x),1);
for ii = 1:length(x)
    nodevec(ii) = node(x(ii));
end
nodevec = nodevec(:);
signnodevec = (-1).^max(m,n)*sign(nodevec);

if ~logFlag
    cvx_quiet(true);
end

if ( m > n)
    for ii = 1:m+1
        TpowersDenom = [TpowersDenom 1./(x-tk(ii))];
    end
    TpowersNum=TpowersDenom;

    % find an initial solution by solving the feasibilty
    % problem
    cvx_begin
    cvx_quiet(true);
    variable a(m+1);
    variable b(m+1);
    subject to
         abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
             signnodevec.*(errPrev*TpowersDenom*b);
         abs(b) <= ones(m+1,1);
         Qmn'*b == 0;
    cvx_end

    aPrev = a; bPrev = b;
    rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
    err = max(abs(fx-rx));

    % perform the DC iteration until convergence
    while abs(abs(err)-abs(errPrev))/abs(err) > 1e-6
        errPrev = err;
        cvx_begin
        variable a(m+1);
        variable b(m+1);
        variable varEps;
        minimize varEps;
        subject to
            abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
                signnodevec.*(err*TpowersDenom*b)+...
                signnodevec.*(varEps*TpowersDenom*bPrev);
            abs(b) <= ones(m+1,1);
            Qmn'*b == 0;
        cvx_end

        aPrev = a; bPrev = b;
        rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
        err = max(abs(fx-rx));
    end

elseif ( m < n )
    for ii = 1:n+1
        TpowersDenom = [TpowersDenom 1./(x-tk(ii))];
    end
    TpowersNum=TpowersDenom;

    % find an initial solution by solving the feasibilty
    % problem
    cvx_begin
    variable a(n+1);
    variable b(n+1);
    subject to
         abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
             signnodevec.*(errPrev*TpowersDenom*b);
         abs(b) <= ones(n+1,1);
         Qmn'*a == 0;
    cvx_end

    aPrev = a; bPrev = b;
    rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
    err = max(abs(fx-rx));

    % perform the DC iteration until convergence
    while abs(abs(err)-abs(errPrev))/abs(err) > 1e-6
        errPrev = err;
        cvx_begin
        variable a(n+1);
        variable b(n+1);
        variable varEps;
        minimize varEps;
        subject to
            abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
                signnodevec.*(err*TpowersDenom*b)+ ...
                signnodevec.*(varEps*TpowersDenom*bPrev);
            abs(b) <= ones(n+1,1);
            Qmn'*a == 0;
        cvx_end

        aPrev = a; bPrev = b;
        rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
        err = max(abs(fx-rx));
    end

else
    for ii = 1:m+1
        TpowersDenom = [TpowersDenom 1./(x-tk(ii))];
    end
    TpowersNum=TpowersDenom;

    % find an initial solution by solving the feasibilty
    % problem
    cvx_begin
    cvx_quiet(true);
    variable a(m+1);
    variable b(m+1);
    subject to
         abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
             signnodevec.*(errPrev*TpowersDenom*b);
         abs(b) <= ones(m+1,1);
    cvx_end

    aPrev = a; bPrev = b;
    rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
    err = max(abs(fx-rx));


    while abs(abs(err)-abs(errPrev))/abs(err) > 1e-6
        errPrev = err;
        cvx_begin
        variable a(m+1);
        variable b(m+1);
        variable varEps;
        minimize varEps;
        subject to
            abs(TpowersNum*a-fx.*(TpowersDenom*b)) <= ...
                signnodevec.*(err*TpowersDenom*b)+ ...
                signnodevec.*(varEps*TpowersDenom*bPrev);
            abs(b) <= ones(m+1,1);
        cvx_end

        aPrev = a; bPrev = b;
        rx = TpowersNum*aPrev./(TpowersDenom*bPrev);
        err = max(abs(fx-rx));
    end
end

[xk,err] = locateExtrema(x,fx-rx);
D = @(zz) 0; N = @(zz) 0;    % form function handle rh = N/D
for ii = 1:length(tk)
   D = @(zz) D(zz) + b(ii)./(zz-tk(ii));
   N = @(zz) N(zz) + a(ii)./(zz-tk(ii));
end
rh = @(zz) reval(zz, tk, N, D, a, b);
if displayFlag
    plot(x,fx-rh(x)),hold on,
    plot(xk,err,'.k','markersize',14),
    title(['Type (', num2str(m), ',', num2str(n), ...
        ') approximation error']), hold off;
    pause(1);
end
err = max(abs(err));
if logFlag
    fprintf('Type (%d,%d) approximation error = %.10f\n', m, n, err);
end

end

function r = reval(zz, xsupport, N, D, wN, wD)
    zv = zz(:);
    r = N(zv)./D(zv);
    ii = find(isnan(r));
    for jj = 1:length(ii)
        if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == xsupport) )
            % r(NaN) = NaN is fine.
            % The second case may happen if r(zv(ii)) = 0/0 at some point.
        else
            % Clean up values NaN = inf/inf at support points.
            % Find the corresponding node and set entry to correct value:
            pos = zv(ii(jj)) == xsupport;
            r(ii(jj)) = -wN(pos)./wD(pos);
        end
    end
    r = reshape(r, size(zz));
end


function [s,es] = locateExtrema(r, er)

% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest
        % value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];
        es = [es ; er(i)];
    end
end

maxVal = max(abs(es));
idx = find(abs(abs(es)-maxVal)/maxVal < 1e-2);
s = s(idx);
es = es(idx);

end


function nxk = refGen(x, xk, n, symType)

xx = linspace(-1,1,length(xk));

if ( symType == 0 )

    nxk = pwiselin(xx, xk, n, linspace(-1,1,n));

% handling of even symmetries
elseif ( symType == 1 )

    halfSize = length(xx)/2;
    halfn = n/2;


    if (xk(1) == x(1))
        nxk = pwiselin(xx(halfSize+1:end), xk(halfSize+1:end),halfn, ...
            linspace(xx(halfSize+1),xx(end),halfn));
        nxk = [nxk; -nxk(2:end); x(1)];
        nxk = sort(nxk,'ascend');
    elseif (xk(end) == x(end))
        nxk = pwiselin(xx(1:halfSize), xk(1:halfSize),halfn, ...
            linspace(xx(1),xx(halfSize),halfn));
        nxk = [nxk; -nxk(1:end-1); x(end)];
        nxk = sort(nxk,'ascend');
    else
        nxk = pwiselin(xx, xk, n, linspace(-1,1,n));
    end

% handling of odd symmetries
else

halfSize = (length(xx)-1) / 2;
halfn = (n-1)/2;

    if ( xk(1) == x(1) )
        nxk = pwiselin(xx(halfSize+2:end), xk(halfSize+2:end),halfn, ...
            linspace(xx(halfSize+2),xx(end),halfn));
        nxk = [nxk; -nxk(1:end); f.domain(1)];
        nxk = sort(nxk,'ascend');
    elseif ( xk(end) == f.domain(end) )
        nxk = pwiselin(xx(1:halfSize+1), xk(1:halfSize+1),halfn, ...
            linspace(xx(1),xx(halfSize+1),halfn));
        nxk = [nxk; -nxk(1:end); x(end)];
        nxk = sort(nxk,'ascend');
    else
        nxk = pwiselin(xx, xk, n, linspace(-1,1,n));
    end

end

end

function yi = pwiselin(xd, yd, ni, xi)

  nd = length(xd);
  xd = xd(:);
  yd = yd(:);
  xi = xi(:);

  if ( nd == 1 )
    yi(1:ni,1) = yd;
    return
  end

  [~, ~,k] = histcounts(xi, xd);

  k ( k == 0 ) = 1;
  k ( k == nd ) = nd - 1;

  t = ( xi - xd(k,1) ) ./ ( xd(k+1,1) - xd(k,1) );
  yi = ( 1 - t ) .* yd(k) + t .* yd(k+1);

  return
end
