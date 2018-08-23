# MAT
math


function [x,fx,xspan]= ShiWeiFa_falsi(fun,a,b,ep)
% FASLI 试位法求根
% 优点计算过程简单收敛可保证 ，对函数性质要求低，只要求连续即可。缺点是收敛速度较慢
% x=b-((b-a)/(f(b)-f(a)))*fb
if nargin==3
    ep=1e-8;   %默认精度
end
fa=fun(a);
fb=fun(b);  % 区间端点处的函数值
xspan=[a,b];
if fa*fb>=0
    error('二分法求根必须要求两端点处函数值.')   %不满足二分法条件
end
while abs(b-a)>ep
    x=b-(b-a)/(fb-fa)*fb;
    fx=fun(x);
    if fx*fa<0;
        b=x;
        fb=fx;
    else 
        a=x;
        fa=fx;
    end
    xspan=[xspan ;a,b];
    if abs(fx)<ep
        break   %若区间端点处的函数值充分小，则认为此端点是方程的根
    end
end

%ep. 利用试位法求函数f(x)=sin(x)-x-x^2+0.5在区间[0,1]内的零点

% f=@(x)sin(x)-x-x.^2+0.5;
% [x,fval,xspan]=ShiWeiFa_falsi(f,0,1,1e-6)

% fplot(f,[0,1])
% hold on
% plot(xlim,[0,0],'--',x,fval,'k*','markersize',10)
% title({['二分法求得的根：x*=',num2str(x)],   ['迭代次数:n=',int2str(size(xspan,1)-1)]   })


2.

function [x,fx,xspan]= ErFenFa_bisect(fun,a,b,ep)
% BISECT 二分法求根
% 优点计算过程简单收敛可保证 ，对函数性质要求低，只要求连续即可。缺点是收敛速度较慢

if nargin==3
    ep=1e-8;   %默认精度
end
fa=fun(a);
fb=fun(b);  % 区间端点处的函数值
xspan=[a,b];
if fa*fb>=0
    error('二分法求根必须要求两端点处函数值.')   %不满足二分法条件
end
while abs(b-a)>ep
    x=(a+b)/2;
    fx=fun(x);
    if fx*fa<0;
        b=x;
        fb=fx;
    else 
        a=x;
        fa=fx;
    end
    xspan=[xspan ;a,b];
    if abs(fx)<ep
        break   %若区间端点处的函数值充分小，则认为此端点是方程的根
    end
end


%ep. 利用试位法求函数f(x)=sin(x)-x-x^2+0.5在区间[0,1]内的零点

% f=@(x)sin(x)-x-x.^2+0.5;
% [x,fval,xspan]=ShiWeiFa_falsi(f,0,1,1e-6)

% fplot(f,[0,1])
% hold on
% plot(xlim,[0,0],'--',x,fval,'k*','markersize',10)
% title({['二分法求得的根：x*=',num2str(x)],   ['迭代次数:n=',int2str(size(xspan,1)-1)]   })

3.

function varargout = fixpt( phifun,x0,ep,maxiter )
%FIXPT 不动点迭代法求方程的根
%   此处显示详细说明
    if nargin<4
        maxiter=500;    % 如果没输最大迭代次数maxiter ，默认为500
    end
    if nargin<3         %如果没输最大允许误差ep ，默认为10^-8
        ep=1e-8;
    end
    iter=1;                 %迭代次数
    xs (iter,1)=x0;         %迭代序列初始值
    exitflag=1;             %迭代发散标志，1表示迭代收敛，0表示迭代发散
    while exitflag 
        x1=phifun(x0);      %迭代计算函数值
        xs(iter+1,1)=x1;    %将迭代值依此存入迭代序列中
        if abs(x1-x0)<=ep   %前后两次迭代值差的绝对值在允许的误差范围内
            break           %跳出循环
        end
        x0=x1;              %更新迭代初始值
        iter=iter+1;
        if iter>maxiter     %若迭代次数大于最大迭代次数，则认为迭代发散，即根不可靠
            exitflag=0;
            break
        end
    end
    [varargout{1:4}]=deal(x1,...  %第1个输出参数为函数零点
    exitflag,...                 %第2个输出参数为迭代收敛标志
    iter,...                     %第3个输出参数为迭代次数
    xs);                         %第4个输出参数为迭代序列
end

%ep. 利用不动点迭代法求解函数f(x)=exp(-x)-x在0.1附近的零点

%首先将方程f(x)=0写成不动点方程的形式： exp(-x)=x ; 然后调用fixpt函数求解即可，编写如下语句：
% phi=@(x)exp(-x);                            %定义不动点函数
% [x,exitflag,iter,Xs]=fixpt(phi,0.1,1e-6);          %不动点迭代法求根
% fplot(@(x)[phi(x),x],[0,1])                 %绘制不动点函数的曲线及y=x的图形
% hold on                                     %图形保持
% Xr=repelem(Xs,2);                           %将Xs的每个元素均复制两次
% for k=1:2*iter 
%     plot(Xr(k:k+1),Xr(k+1:k+2),'k')         %动态绘制每一段线段       
%     pause(0.1)                              %暂停0.1秒
% end
% title({['不动点迭代法求得的根：x*=',num2str(x)],['迭代次数：n=',int2str(iter)]})

4.

function varargout=newtonDieDai(fun,x0,ep,maxiter)
% NEWTOM  牛顿法求方程的根
if nargin<4
    maxiter=500;   %如果没输最大迭代次数maxiter ，默认为500
end
if nargin<3
    ep=1e-8;      %默认允许误差
end
if ~isscalar(fun)
    dfun=fun{2};                          %导函数匿名函数形式
    fun=fun{1};                           %函数的匿名函数形式
else
    if isa(fun,'sym')                     %函数以符号表达式的形式给出
        dfun=matlabFunction(diff(fun));   %导函数匿名函数形式
        fun=matlabFunction(fun);          %函数的匿名函数形式
    elseif isa(fun,'function_handle')     %函数以匿名函数或函数句柄的形式给出
        dfun=matlabFunction(diff(sym(@(x)fun(x)))); %导函数匿名函数形式
    end
end
iter=1;        %迭代次数
xs(iter,1)=x0; %迭代序列初始值
exitflag=1;    % %迭代发散标志，1表示迭代收敛，0表示迭代发散
x1=nan;        %预置x1的初值
while exitflag
    fx=fun(x0);     %计算x0处的函数值
    dfx=dfun(x0);   %计算x0处的导数值
    if abs(dfx)<=eps|| iter>maxiter  %若导数为0或迭代次数大于最大迭代次数
        exitflag=0;    % 认为迭代发散，即根不可靠
        break         %退出循环
    end
    x1=x0-fx/dfx;     %牛顿迭代计算
    xs(iter+1,1)=x1;  %将迭代值依此存入迭代序列中
    if abs(x1-x0)<=ep %前后两次迭代值差的绝对值在允许的误差范围内
        break         %跳出循环
    end
    x0=x1;
    iter=iter+1;
end
[varargout{1:5}]=deal(x1,... %第1个输出参数为函数零点
       fun(x1),...           %第2个输出参数为函数零点的函数值
       exitflag,...          %第3个输出参数为迭代收敛标志
       iter,...              %第4个输出参数为迭代次数
       xs);                  %第5个输出参数为迭代序列

% 利用牛顿法求解函数f(x)=exp(x)-x-5在3.8附近的零点

% fun=@(x)exp(x)-x-5;
% [x,fval,exitflag,iter,Xs]=newtonDieDai(fun,3.8,1e-6)
% fplot(@(x)[fun(x),zeros(size(x))],[1,4])
% hold on 
% Xr=repelem(Xs,2);
% Yr=reshape([zeros(size(Xs)),fun(Xs)].',[ ],1);  %在fun(Xs)的各元素前面插入0
% for k=1:2*iter+1
%     plot(Xr(k:k+1),Yr(k:k+1),'k')  %动态绘制每一段线段
%     pause (0.2)                    %暂停0.2秒
% end
% plot(x,fval,'k*','markersize',6)   %绘制零点
% title({['牛顿法求得的根：x*=',num2str(x)],...
%     ['迭代次数：n=',int2str(iter)]})

%给定一个复数方程f(z)=0,在z平面内任取一点作为初值，利用牛顿法求f(z)=0的根，要求对使迭代收敛的初始点进行着色,直至
%z平面内的所有点都迭代完毕为止
fun=@(z)z.^3-1;      %复数函数
dfun=@(z)3*z.^2;     %导数
xspan=[-2,2]; yspan=[-2,2]; %绘图区域范围
x=linspace(xspan(1),xspan(end),512); %x方向等间隔采样点
y=linspace(yspan(1),xspan(end),512); %x方向等间隔采样点
[X,Y]=meshgrid(x,y);                 %利用向量生成矩阵
z=X+1i*Y;                            %创建复数矩阵
W=zeros(size(z));                    
for r=1:size(z,1)
    for c=1:size(z,2)
        [~,~,exitflag,iter]=newtonDieDai({fun,dfun},z(r,c),1e-5,20);
        if exitflag           %迭代法求复根，若迭代收敛
            W(r,c)=iter;      %将收敛点的迭代次数存入矩阵W中
        end
    end
end
image(x,y,W)     %绘制图像
colormap hsv   %设置色图矩阵

5.
function y = linspace(d1, d2, n)
%LINSPACE Linearly spaced vector.
%   LINSPACE(X1, X2) generates a row vector of 100 linearly
%   equally spaced points between X1 and X2.
%
%   LINSPACE(X1, X2, N) generates N points between X1 and X2.
%   For N = 1, LINSPACE returns X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   See also LOGSPACE, COLON.

%   Copyright 1984-2016 The MathWorks, Inc.

if nargin == 2
    n = 100;
else
    n = floor(double(n));
end
if ~isscalar(d1) || ~isscalar(d2) || ~isscalar(n)
    error(message('MATLAB:linspace:scalarInputs'));
end
n1 = n-1;
c = (d2 - d1).*(n1-1); %check intermediate value for appropriate treatment
if isinf(c)
    if isinf(d2 - d1) %opposite signs overflow
        y = d1 + (d2./n1).*(0:n1) - (d1./n1).*(0:n1);
    else
        y = d1 + (0:n1).*((d2 - d1)./n1);
    end
else
    y = d1 + (0:n1).*(d2 - d1)./n1;
end
if ~isempty(y)
    if d1 == d2
        y(:) = d1;
    else
        y(1) = d1;
        y(end) = d2;
    end
end

6.
function c = nchoosek(v,k)
%NCHOOSEK Binomial coefficient or all combinations.
%   NCHOOSEK(N,K) where N and K are non-negative integers returns N!/K!(N-K)!.
%   This is the number of combinations of N things taken K at a time.
%   When a coefficient is large, a warning will be produced indicating
%   possible inexact results. In such cases, the result is only accurate
%   to 15 digits for double-precision inputs, or 8 digits for single-precision
%   inputs.
%
%   NCHOOSEK(V,K) where V is a vector of length N, produces a matrix
%   with N!/K!(N-K)! rows and K columns. Each row of the result has K of
%   the elements in the vector V. This syntax is only practical for
%   situations where N is less than about 15.
%
%   Class support for inputs N,K:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   Class support for inputs V:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%      logical, char
%
%   See also PERMS.

%   Copyright 1984-2013 The MathWorks, Inc.

if ~isscalar(k) || k < 0 || ~isreal(k) || k ~= round(k)
    error(message('MATLAB:nchoosek:InvalidArg2'));
end

if ~isvector(v)
    error(message('MATLAB:nchoosek:InvalidArg1'));
end

% the first argument is a scalar integer
if isscalar(v) && isnumeric(v) && isreal(v) && v==round(v) && v >= 0
    % if the first argument is a scalar, then, we only return the number of
    % combinations. Not the actual combinations.
    % We use the Pascal triangle method. No overflow involved. c will be
    % the biggest number computed in the entire routine.
    %
    n = v;  % rename v to be n. the algorithm is more readable this way.
    if isinteger(n)
        if ~(strcmp(class(n),class(k)) || isa(k,'double'))
            error(message('MATLAB:nchoosek:mixedIntegerTypes'))
        end
        classOut = class(n);
        inttype = true;
        int64type = isa(n,'int64') || isa(n,'uint64');
    elseif isinteger(k)
        if ~isa(n,'double')
            error(message('MATLAB:nchoosek:mixedIntegerTypes'))
        end
        classOut = class(k);
        inttype = true;
        int64type = isa(k,'int64') || isa(k,'uint64');
    else % floating point types
        classOut = superiorfloat(n,k);
        inttype = false;
        int64type = false;
    end
    
    if k > n
        error(message('MATLAB:nchoosek:KOutOfRange'));
    elseif ~int64type && n > flintmax
        error(message('MATLAB:nchoosek:NOutOfRange'));
    end
    
    if k > n/2   % use smaller k if available
        k = n-k;
    end
    
    if k <= 1
        c = n^k;
    else
        if int64type
            % For 64-bit integers, use an algorithm that avoids
            % converting to doubles
            c = binCoef(n,k,classOut);
        else
            % Do the computation in doubles.
            nd = double(n);
            kd = double(k);
            
            nums = (nd-kd+1):nd;
            dens = 1:kd;
            nums = nums./dens;
            c = round(prod(nums));
            
            if ~inttype && c > flintmax(classOut)
                warning(message('MATLAB:nchoosek:LargeCoefficient', ...
                    sprintf( '%e', flintmax(classOut) ), floor(log10(flintmax(classOut)))));
            end
            % Convert answer back to the correct type
            c = cast(c,classOut);
        end
    end
    
else
    % the first argument is a vector, generate actual combinations.
    
    n = length(v);
    if iscolumn(v)
        v = v.';
    end
    
    if n == k
        c = v;
    elseif n == k + 1
        c = repmat(v,n,1);
        c(1:n+1:n*n) = [];
        c = reshape(c,n,k);
    elseif k == 1
        c = v.';
    elseif k == 0
        c = zeros(1,0,class(v));
    elseif n < 17 && (k > 3 || n-k < 4)
        tmp = uint16(2^n-1):-1:2;
        x = bsxfun(@bitget,tmp.',n:-1:1);
        
        idx = x(sum(x,2) == k,:);
        nrows = size(idx,1);
        [rows,~] = find(idx');
        c = reshape(v(rows),k,nrows).';
    else
        [~,maxsize] = computer;
        % Error if output dimensions are too large
        if k*nchoosek(n,k) > maxsize
            error(message('MATLAB:pmaxsize'))
        end
        c = combs(v,k);
    end
end

end

%----------------------------------------
function c = binCoef(n,k,classOut)
% For integers, compute N!/((N-K)! K!) using prime factor cancellations

numerator = cast((n-k+1):n,classOut);
for denominator = k:-1:1
    F = factor(denominator);
    for whichfactor = 1:numel(F)
        thefactor = F(whichfactor);
        largestMultiple = find(mod(numerator,thefactor) == 0, 1, 'last');
        numerator(largestMultiple) = numerator(largestMultiple)/thefactor;
    end
end
c = prod(numerator,'native');
end

%----------------------------------------
function P = combs(v,m)
%COMBS  All possible combinations.
%   COMBS(1:N,M) or COMBS(V,M) where V is a row vector of length N,
%   creates a matrix with N!/((N-M)! M!) rows and M columns containing
%   all possible combinations of N elements taken M at a time.
%
%   This function is only practical for situations where M is less
%   than about 15.

v = v(:).'; % Make sure v is a row vector.
n = length(v);
if n == m
    P = v;
elseif m == 1
    P = v.';
else
    P = [];
    if m < n && m > 1
        for k = 1:n-m+1
            Q = combs(v(k+1:n),m-1);
            P = [P; [v(ones(size(Q,1),1),k) Q]]; %#ok
        end
    end
end
end


function [a,b] = polyder(u,v)
%POLYDER Differentiate polynomial.
%   POLYDER(P) returns the derivative of the polynomial whose
%   coefficients are the elements of vector P.
%
%   POLYDER(A,B) returns the derivative of polynomial A*B.
%
%   [Q,D] = POLYDER(B,A) returns the derivative of the
%   polynomial ratio B/A, represented as Q/D.
%
%   Class support for inputs u, v:
%      float: double, single
%
%   See also POLYINT, CONV, DECONV.

%   Copyright 1984-2007 The MathWorks, Inc.

if nargin < 2, v = 1; end

u = u(:).'; v = v(:).';
nu = length(u); nv = length(v);
if nu < 2, up = 0; else up = u(1:nu-1) .* (nu-1:-1:1); end
if nv < 2, vp = 0; else vp = v(1:nv-1) .* (nv-1:-1:1); end
a1 = conv(up,v); a2 = conv(u,vp);
i = length(a1); j = length(a2); z = zeros(1,abs(i-j));
if i > j, a2 = [z a2]; elseif i < j, a1 = [z a1]; end
if nargout < 2, a = a1 + a2; else a = a1 - a2; end
f = find(a ~= 0);
if ~isempty(f), a = a(f(1):end); else a = zeros(superiorfloat(u,v)); end
b = conv(v,v);
f = find(b ~= 0);
if ~isempty(f), b = b(f(1):end); else b = zeros(class(v)); end
%  The vector may be too long when the polynomial coefficients
%  include NaN, so trim the leading element.
if length(a) > max(nu + nv - 2,1)
  a = a(2:end);
end
