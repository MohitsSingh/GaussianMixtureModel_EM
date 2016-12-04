
clear all;

%load "digit.mat";


[dx,n] = settings()
m = 3;

tws = zeros(m,1);
remn = n;
for j=1:m-1
  tws(j) = randi(round(remn*0.75));
  remn = remn - tws(j);
end
tws(m) = remn;

tmus = randn(m,dx)*10;%[-5,0,8]';%
tsigmas = zeros(dx,dx,m); %rand(d,d,m);
for k = 1:m
  tsigmas(:,:,k) = speye(dx,dx);
end

xx = sqrtm(tsigmas(:,:,1)) * randn(dx,tws(1)) + tmus(1,:)';
for k = 2:m
  xx = horzcat(xx,sqrtm(tsigmas(:,:,k))*randn(dx,tws(k)) + tmus(k,:)');
end
xx = xx';
tws = tws ./ n;

h3 = figure;
histogram(xx(:,1),100);

%histogram(xx(:,1));
x = -10:0.1:10;
dx = size(x);
dx = dx(2);
y1 = zeros(1,dx);
y2 = zeros(1,dx);
y3 = zeros(1,dx);
y4 = zeros(1,dx);

size(tws)
size(tmus)
size(xx)
size(tsigmas)
for xi = 1:dx
  y1(1,xi) = tws(1)*N2(x(xi),tmus(1,:),tsigmas(:,:,1));
  y2(1,xi) = tws(2)*N2(x(xi),tmus(2,:),tsigmas(:,:,2));
  y3(1,xi) = tws(3)*N2(x(xi),tmus(3,:),tsigmas(:,:,3));
  y4(1,xi) = y1(1,xi) + y2(1,xi) +y3(1,xi);
end
h1=figure;
fig = plot(x,y1,x,y2,x,y3,x,y4);
title('å≥ÉfÅ[É^')

pause;

%return

[ws,mus,sigmas] = MixMLE(xx,3,1.0e-8);


yy1 = zeros(1,dx);
yy2 = zeros(1,dx);
yy3 = zeros(1,dx);
yy4 = zeros(1,dx);
for xi = 1:dx
  yy1(1,xi) = ws(1)*N2(x(xi),mus(1,:),sigmas(:,:,1));
  yy2(1,xi) = ws(2)*N2(x(xi),mus(2,:),sigmas(:,:,2));
  yy3(1,xi) = ws(3)*N2(x(xi),mus(3,:),sigmas(:,:,3));
  yy4(1,xi) = yy1(1,xi) + yy2(1,xi) +yy3(1,xi);
end
h2=figure;
fig2 = plot(x,yy1,x,yy2,x,yy3,x,yy4);
title('ìöÇ¶')
pause;
close(h1);
close(h2);

%normpdf
function y = N2(ax,mu,sigma)
  d = size(ax);
  mu = mu(1);
  sigma = sigma(1);
  d = d(1);
  y= 1 / ( (2 * pi)^(d/2) * det(sigma)^(1/2)) * exp(-1/2 * (ax-mu)' * pinv(sigma) * (ax-mu));
end