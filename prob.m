
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
      dmus = dmus' * dmus;
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
subplot(1,3,1);
hist1 = histogram(xx(1,:),100);
title('�W�{�̕��z')

ax = gca;
xlim = ax.XLim;
x = xlim(1):0.2:xlim(2);
%�`��f�[�^����
y1(1,:) = (tws(1)*mvnpdf(x',tmus(:,1)',tsigmas(:,:,1)))';
y2(1,:) = (tws(2)*mvnpdf(x',tmus(:,2)',tsigmas(:,:,2)))';
y3(1,:) = (tws(3)*mvnpdf(x',tmus(:,3)',tsigmas(:,:,3)))';
y4(1,:) = y1(1,:) + y2(1,:) + y3(1,:);

subplot(1,3,2);
fig = plot(x,y1,x,y2,x,y3,x,y4);
title('�^�̕��z')

pause;

%===================
%�Ŗސ���̎��s
%===================

%���Z�ߒ��̃O���t�������ꍇ��MixMLE
[ws,mus,sigmas] = MixMLE(xx,3,1.0e-4);

%===================
%�Ŗސ���ʂŃO���t�쐬
%===================
yy1(1,:) = (ws(1)*mvnpdf(x',mus(:,1)',sigmas(:,:,1)))';
yy2(1,:) = (ws(2)*mvnpdf(x',mus(:,2)',sigmas(:,:,2)))';
yy3(1,:) = (ws(3)*mvnpdf(x',mus(:,3)',sigmas(:,:,3)))';
yy4(1,:) = yy1(1,:) + yy2(1,:) + yy3(1,:);
gcf = h1;
subplot(1,3,3);
fig2 = plot(x,yy1,x,yy2,x,yy3,x,yy4);
title('�Ŗސ��茋��')

pause;
close(h1);
