%ガウス混合モデルに対して最尤推定を行い最適なパラメータを返す関数
%引数
% xs d*n行列 各行にd次元標本が入っている
% m 混合数
% e 収束判定に使う正数 etaの差がこれより小さくなったら終了する
%返り値
% ws = m * 1 の列ベクトル
% mus = d*m の行列 muの行ベクトルが並んでいる
% sigmas = d*d*n 分散共分散行列がn枚並んでいる
%グラフの描画あり
function [ws,mus,sigmas] = MixMLE(xs,m,e)
  [d,n] = size(xs);
  %初期値設定
  ws = zeros(m,1);
  gamma = rand(m,1);
  gsum  = 0;
  for j=1:m
    ws(j) = exp(gamma(j));
    gsum = gsum + ws(j);
  end
  ws = ws ./ gsum;
  mus = rand(d,m);
  sigmas = rand(d,d,m);
  eta = zeros(n,m);
  
  drawfreq = 5;
  count = drawfreq;
  h1 = figure;
  while 1
    count = count + 1;
    lasteta = eta;
    %Eステップ
    for i=1:n
      x = xs(:,i);
      eta(i,:) = (ws .* mvnpdf(repmat(x,[1,m])',mus',sigmas))';
      eta(i,:) = eta(i,:) ./ sum(eta(i,:),2);
    end
    
    %収束判定
    diff = sum(sum(~((eta - lasteta) < e)))
    if diff == 0
      break;
    end
    
    %drawfreq回に１回グラフを描画
    if  drawfreq < count
      count = 0;
      x = -10:0.5:10;
      y1(1,:) = (ws(1)*mvnpdf(x',mus(:,1)',sigmas(:,:,1)))';
      y2(1,:) = (ws(2)*mvnpdf(x',mus(:,2)',sigmas(:,:,2)))';
      y3(1,:) = (ws(3)*mvnpdf(x',mus(:,3)',sigmas(:,:,3)))';
      y4(1,:) = y1(1,:) + y2(1,:) + y3(1,:);
      plot(x,y1,x,y2,x,y3,x,y4);
      title('計算過程')
      drawnow
    end
    
    %Mステップ
    temp = (sum(eta))';
    %重みの計算
    ws = temp ./ n;
    %平均の計算
    temp1 = repmat(eta,[1,1,d]);%n,m,d
    temp1 = permute(temp1,[3,2,1]);%
    temp2 = repmat(xs,[1,1,m]);%d,n,m
    temp2 = permute(temp2,[1,3,2]);%d,m,n
    mus = sum((temp1.*temp2),3);%d,m
    mus = mus ./ (repmat(temp,[1,d]))';
    %分散共分散行列の計算
    temp3 = repmat(xs,[1,1,m]);%d*n*m
    temp3 = permute(temp3,[1,3,2]);%dmn
    temp4 = repmat(mus,[1,1,n]);%d*m*n
    temp5 = temp3 - temp4;%dmn
    temp6 = repmat(temp5,[1,1,1,d]);%dmnd
    temp6 = permute(temp6,[1,4,2,3]);%ddmn
    temp7 = permute(temp6,[2,1,3,4]);%ddmn
    temp8 = temp6 .* temp7;%d*d*m*n
    temp10 = repmat(eta,[1,1,d,d]);%n*m*d*d
    temp10 = permute(temp10,[3,4,2,1]);%d*d*m*n
    sigmas = sum((temp8.*temp10),4);%d*d*m
    temp11 = (repmat(temp,[1,d,d]));%m*d*d
    temp11 = permute(temp11,[2,3,1]);%d*d*m
    sigmas = sigmas ./ temp11;%d*d*m
  end
  close(h1);
end