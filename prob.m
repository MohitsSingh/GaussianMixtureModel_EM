
clear

%�ݒ�ǂݍ���
[d,n] = settings()
m = 3;

%===================
%�^�̃f�[�^����
%===================

%�^�̏d��
tws = zeros(m,1);
remn = n;
for j=1:m-1
  tws(j) = randi(round(remn*0.75));
  remn = remn - tws(j);
end
tws(m) = remn;

%�^�̕���
%tmus = rand(d,m)*20 - 10;%[-5,0,8];%[-5,0,8]';%
%tmus�������悤�Ɍv�Z����
tmus(:,1) = rand(d,1)*20-10;
for j = 2:m
  while 1
    tmus(:,j) = rand(d,1)*20-10;
    done = 1;
    for i = 1:j-1
      dmus = tmus(:,i) - tmus(:,j);
      dmus = dmus' * dmus
      if dmus < 4
        done = 0;
        break
      end
    end
    if done == 1
      break;
    end
  end
end

%�^�̕��U�����U�s��
tsigmas = zeros(d,d,m); %rand(d,d,m);
for k = 1:m
  tsigmas(:,:,k) = speye(d,d);
end

%�W�{����
xx = sqrtm(tsigmas(:,:,1)) * randn(d,tws(1)) + tmus(:,1);
for k = 2:m
  xx = horzcat(xx,sqrtm(tsigmas(:,:,k))*randn(d,tws(k)) + tmus(:,k));
end
tws = tws ./ n;

%===================
%�O���t�`��
%===================

h1 = figure;
h1.Position = [10 309 1300 400];
gcf = h1;
subplot(1,3,1);
hist1 = histogram(xx(1,:),100);
title('�W�{�̕��z')

ax = gca;
xlim = ax.XLim;
x = xlim(1):0.2:xlim(2);
xc = size(x);
xc = xc(2);
y1 = zeros(1,xc);
y2 = zeros(1,xc);
y3 = zeros(1,xc);
y4 = zeros(1,xc);
for xi = 1:xc
  y1(1,xi) = tws(1)*N2(x(:,xi),tmus(:,1),tsigmas(:,:,1));
  y2(1,xi) = tws(2)*N2(x(:,xi),tmus(:,2),tsigmas(:,:,2));
  y3(1,xi) = tws(3)*N2(x(:,xi),tmus(:,3),tsigmas(:,:,3));
  y4(1,xi) = y1(1,xi) + y2(1,xi) +y3(1,xi);
end
subplot(1,3,2);
fig = plot(x,y1,x,y2,x,y3,x,y4);
title('�^�̕��z')

pause;

%===================
%�Ŗސ���̎��s
%===================

%���Z�ߒ��̃O���t�������ꍇ��MixMLE
[ws,mus,sigmas] = MixMLE1(xx,3,1.0e-4);

%===================
%�Ŗސ���ʂŃO���t�쐬
%===================

yy1 = zeros(1,xc);
yy2 = zeros(1,xc);
yy3 = zeros(1,xc);
yy4 = zeros(1,xc);
for xi = 1:xc
  yy1(1,xi) = ws(1)*N2(x(:,xi),mus(:,1),sigmas(:,:,1));
  yy2(1,xi) = ws(2)*N2(x(:,xi),mus(:,2),sigmas(:,:,2));
  yy3(1,xi) = ws(3)*N2(x(:,xi),mus(:,3),sigmas(:,:,3));
  yy4(1,xi) = yy1(1,xi) + yy2(1,xi) +yy3(1,xi);
end
gcf = h1;
subplot(1,3,3);
fig2 = plot(x,yy1,x,yy2,x,yy3,x,yy4);
title('�Ŗސ��茋��')

pause;
close(h1);




%���K���z�̊m����Ԃ��֐�
%ax = d*1
%mu = d*1
%sigma = d*d
function y = N2(ax,mu,sigma)
  d = size(ax);
  d = d(1);
  y= 1 / ( (2 * pi)^(d/2) * det(sigma)^(1/2)) * exp(-1/2 * (ax-mu)' * pinv(sigma) * (ax-mu));
end