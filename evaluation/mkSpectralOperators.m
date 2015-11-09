function varargout = mkSpectralOperators(varargin)

%computes [D_i,D^2_i,x_i,Int_i].  The D's and Int's are matricies formed from
%kronecker products.  The input varargin{i} = {M_i,'kind'} where M_i number of
%points in ith dirction and kind = cheb or kind = fourier.  Cheby grids are assumed to be x = cos(pi*(
%2*n+1)/(2*M) ) and fourier grids are -pi + 2*pi*n/N where N is odd

%constants of integration for cheb grids are fixed at x(end) whereas for
%fourier grids they are fixed at x(1)

Ndims = length(varargin);

d = cell(1,Ndims);
d2 = cell(1,Ndims);
X = cell(1,Ndims);

for n = 1:Ndims
    M = varargin{n}{1};
    kind = varargin{n}{2};
    if strcmp(kind,'cheb')
        [dd,dd2,x] = cheb(M);
    elseif strcmp(kind,'cheb2')
        [dd,dd2,x] = cheb2(M);
    elseif strcmp(kind,'fourier')
        if floor( (M-1)/2 ) ~= (M-1)/2
            error('must use odd number of points for fourier')
        end
        [dd,dd2,x] = fourier((M-1)/2);
    else
        error('unknown kind')
    end
    d{n} = dd;
    d2{n} = dd2; 
    X{n} = x;  
end
    

for n = 1:Ndims
    
    D = 1;
    D2 = 1;
    Int = 1;
    
    for m = Ndims:-1:1
        
        M = varargin{m}{1}; %size of rightmost array
        kind = varargin{m}{2};
        
        if m == n %index equals grad number
            D = kron(D,d{m});
            D2 = kron(D2,d2{m});
            
            if strcmp(kind,'cheb')  
                L = d{m}; L(M,:) = zeros(M,1); L(M,M) = 1; %sets constant of integration to Mth value of function
                Res = zeros(M,M); Res(:,M) = ones(1,M);
                Int = sparse(kron(Int,L\eye(M) - Res));
            elseif strcmp(kind,'cheb2')
                L = d{m}; L(M,:) = zeros(M,1); L(M,M) = 1; %sets constant of integration to Mth value of function
                Res = zeros(M,M); Res(:,M) = ones(1,M);
                Int = sparse(kron(Int,L\eye(M) - Res));
            elseif strcmp(kind,'fourier') 
                L = d{m}; L(1,:) = zeros(M,1); L(1,1) = 1; %sets constant of integration to 1st value of function
                Res = zeros(M,M); Res(:,1) = ones(1,M);
                Int = sparse(kron(Int,L\eye(M) - Res));
            end
        else
            D = kron(D,speye(M));
            D2 = kron(D2,speye(M));
            Int = kron(Int,speye(M));    
        end
        
    end
    
    varargout{n} = D;
    varargout{Ndims+n} = D2;
    varargout{2*Ndims + n} = X{n};
    varargout{3*Ndims + n} = Int;   
end

   

function [D,D2,x] = cheb(N)

k = 1:N; k = k(:);

x = cos( pi*(2*k - 1)/(2*N) );

[n,m] = ndgrid(k,k);

%D = (2*cos((m + n)*pi).*sin((pi - 2*m*pi)/(2*N)))./(-sin(((m - n)*pi)/N) + sin(((-1 + m + n)*pi)/N) + sin((pi - 2*n*pi)/N));

D = (2*(-1).^(n+m).*sin((pi - 2*m*pi)/(2*N)))./(-sin(((m - n)*pi)/N) + sin(((-1 + m + n)*pi)/N) + sin((pi - 2*n*pi)/N));

I = logical(eye(N)); D(I) = 0;

%D(I) = (cos((pi - 2*m(I)*pi)/(2*N)).*csc(((-1 + 2*m(I))*pi)/(2*N)).^2)/2;

D = D - diag(sum(D,2));

D2 = -((-1).^(m + n).*(1 + cos(((m - n).*pi)/N) + cos(((-1 + m + n).*pi)/N) - 3.*cos((pi - 2.*n.*pi)/N)).*csc(((m - n).*pi)/(2.*N)).^2.*csc(((-1 + 2.*n).*pi)/(2.*N)).^3.*sec(((-1 + m + n + N).*pi)./(2.*N)).^2.*sin(((-1 + 2.*m).*pi)/(2.*N)))/8;

D2(I) = 0;

D2 = D2 - diag(sum(D2,2));

function [D,D2,x] = fourier(M)

k = 0:2*M; k = k(:);

x = -pi + 2*pi*k/(2*M+1);

[n,m] = ndgrid(k,k);

D = ((-1).^(m + n).*csc(((-m + n)*pi)/(1 + 2*M)))/2;

I = logical(eye(2*M+1)); 

D(I) = 0;

D2 = (csc(((m - n).*pi)./(1 + 2.*M)).^3.*((1 + M).^2.*sin(((-1 + 2.*M).*(m - n).*pi)./(1 + 2.*M)) + M.^2.*sin(((3 + 2.*M).*(m - n).*pi)./(1 + 2.*M))))./(4 + 8.*M);

D2(I) = -(M*(1 + M))/3;

function [D,D2,x] = cheb2(N)

if N==0, D=0; x=1; return, end
N = N-1;
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));

D2 = zeros(N+1,N+1);

%%% i,j elements
%%%i+j elements
ind = repmat(2:N,N-1,1)+repmat((2:N)',1,N-1);
xx = repmat(x,1,N+1);
D2(2:N,2:N) = 1/2 * (-1).^ind .* (-4 + 2*xx(2:N,2:N).*xx(2:N,2:N) + 2.*xx(2:N,2:N).*xx(2:N,2:N)' ) ./ (1 - xx(2:N,2:N).*xx(2:N,2:N)) ./ (xx(2:N,2:N) - xx(2:N,2:N)').^2 ;
%%i=j elements
diage = logical(eye(N+1,N+1)); diage(1,1) = 0; diage(N+1,N+1) = 0;
xsq = xx.*xx';
D2(diage) = -(5 + N^2 + (1-N^2)*(2.*xsq(diage) - 1))./6./(1-xsq(diage)).^2;

%%other elements:
D2(2:N,1) = (-1).^(2:N)'/2 .* (2 + x(2:N))./(1-x(2:N).^2)./(1-x(2:N));
D2(1,2:N) = ((-1).^(2:N)'/6 .* (5 - 2*N^2 + (1 + 2*N^2).*x(2:N))).*4./((1-x(2:N)).^2);
D2(2:N,N+1) = -(-1).^((2:N)'+N)./2 .* (-2 + x(2:N)).*(1-x(2:N))./(1-x(2:N).^2).^2;
D2(N+1,2:N) = -2.*(-1).^((2:N)'+N)./3 .* (-5 +2*N^2 + (1+2*N^2).*x(2:N)) ./ (1 + x(2:N)).^2;

D2(1,1) = (N^4-1)/15; D2(N+1,N+1) = D2(1,1);
D2(1,N+1) = (-1)^N*(N^2-1)/3; D2(N+1,1) = D2(1,N+1);

%%% comments:
%%%Note minus signs depend on indexing, here indexing in Matlab shifts from 0:N to 1:N+1