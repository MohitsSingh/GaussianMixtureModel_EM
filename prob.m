
clear

%設定読み込み
[d,n] = settings()
m = 3;

%===================
%真のデータ生成
%===================

%真の重み
tws = zeros(m,1);
remn = n;
for j=1:m-1
  tws(j) = randi(round(remn*0.75));
  remn = remn - tws(j);
end
tws(m) = remn;

%真の平均
%tmus = rand(d,m)*20 - 10;%[-5,0,8];%[-5,0,8]';%
%tmusが離れるように計算する
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

%真の分散共分散行列
tsigmas = zeros(d,d,m); %rand(d,d,m);
for k = 1:m
  tsigmas(:,:,k) = speye(d,d);
end

%標本生成
xx = sqrtm(tsigmas(:,:,1)) * randn(d,tws(1)) + tmus(:,1);
for k = 2:m
  xx = horzcat(xx,sqrtm(tsigmas(:,:,k))*randn(d,tws(k)) + tmus(:,k));
end
tws = tws ./ n;

%===================
%グラフ描画
%===================

h1 = figure;
h1.Position = [10 309 1300 400];
subplot(1,3,1);
hist1 = histogram(xx(1,:),100);
title('標本の分布')

ax = gca;
xlim = ax.XLim;
x = xlim(1):0.2:xlim(2);
%描画データ生成
y1(1,:) = (tws(1)*mvnpdf(x',tmus(:,1)',tsigmas(:,:,1)))';
y2(1,:) = (tws(2)*mvnpdf(x',tmus(:,2)',tsigmas(:,:,2)))';
y3(1,:) = (tws(3)*mvnpdf(x',tmus(:,3)',tsigmas(:,:,3)))';
y4(1,:) = y1(1,:) + y2(1,:) + y3(1,:);

subplot(1,3,2);
fig = plot(x,y1,x,y2,x,y3,x,y4);
title('真の分布')

pause;

%===================
%最尤推定の実行
%===================

%演算過程のグラフを示す場合はMixMLE
[ws,mus,sigmas] = MixMLE(xx,3,1.0e-4);

%===================
%最尤推定量でグラフ作成
%===================
yy1(1,:) = (ws(1)*mvnpdf(x',mus(:,1)',sigmas(:,:,1)))';
yy2(1,:) = (ws(2)*mvnpdf(x',mus(:,2)',sigmas(:,:,2)))';
yy3(1,:) = (ws(3)*mvnpdf(x',mus(:,3)',sigmas(:,:,3)))';
yy4(1,:) = yy1(1,:) + yy2(1,:) + yy3(1,:);
gcf = h1;
subplot(1,3,3);
fig2 = plot(x,yy1,x,yy2,x,yy3,x,yy4);
title('最尤推定結果')

pause;
close(h1);
