
%混合ガウスモデルに対して最尤推定を行う
%与えられた標本、混合数に対して最適なパラメータを返す
%xs n*d行列 各行にd次元標本が入っている
%m 混合数
%e 収束判定に使う正数 etaの差がこれより小さくなったら終了する
%返り値
%ws = m * 1 の列ベクトル
%mus = m*d の行列 muの行ベクトルが並んでいる
%sigmas = d*d*n covがn枚並んでいる
function [ws,mus,sigmas] = MixMLE(xs,m,e)
  [n,d] = size(xs);
  %初期値設定
  ws = zeros(m,1);
  gamma = rand(m,1);
  gsum  = 0;
  for j=1:m
    ws(j) = exp(gamma(j));
    gsum = gsum + ws(j);
  end
  ws = ws ./ gsum;
  mus = rand(m,d);
  sigmas = rand(d,d,m);
  %pause;
  
  eta = zeros(n,m);
  lasteta = zeros(n,m);
  
  count = 0;
  lastL = 0;
  while 1
    lasteta = eta;
    %Eステップ
    for i=1:n
      x = xs(i,:)';
      %temp = q(x,ws,mus,sigmas);
      temp = 0;
      for j=1:m
        eta(i,j) = ws(j) * N(x,mus(j,:)',sigmas(:,:,j));
        temp = temp + eta(i,j);
        %eta(i,j) = ws(j) * N(x,mus(j,:)',sigmas(j,:,:)) / temp;
      end
      eta(i,:) = eta(i,:) ./ temp;
    end
    
    %差分計算
    deta = lasteta - eta;
    diff = 0;
    for i=1:n
      for j=1:m
        diff = diff + abs(deta(i,j));
      end
    end
    %deta = deta .* deta;
    %diff = sum(sum(deta));
    diff
    if diff < e
      return
    end
    %diff = sum(sum(~((eta - lasteta) < 0.001)))
    %if diff == 0
    %  return
    %end
    
    %Mステップ
    %temp = m * 1 行列
    temp = (sum(eta))';
    
    %重みの計算
    ws = temp ./ n;
    
    %平均の計算
    temp1 = repmat(eta,[1,1,d]);
    temp1 = permute(temp1,[2,3,1]);
    temp2 = repmat(xs,[1,1,m]);
    temp2 = permute(temp2,[3,2,1]);
    
    mus = sum((temp1.*temp2),3);
    mus = mus ./ (repmat(temp,[1,d]));
    
    
    %Compute Sigma
    temp3 = repmat(xs,[1,1,m]);
    temp3 = permute(temp3,[1,3,2]);
    temp4 = repmat(mus,[1,1,n]);
    temp4 = permute(temp4,[3,1,2]);
    
    %y = n*m*d 行列
    temp5 = temp3 - temp4;
    
    temp6 = repmat(temp5,[1,1,1,d]);
    temp7 = permute(temp6,[1,2,4,3]);
    
    %Y = n*m*d*d行列
    %This is Y n*m*d*d
    temp8 = temp6 .* temp7;
    temp10 = repmat(eta,[1,1,d,d]);
    
    %simgas = m * d* d 
    sigmas = squeeze(sum((temp8.*temp10),1));
    temp11 = (repmat(temp,[1,d,d]));
    sigmas = sigmas ./ temp11;
    
    %sigmas <- d * d * n
    sigmas = permute(sigmas,[2,3,1]);
    
    sigmas
    pause;
    for j=1:m
    %  temp = sum(eta(:,j));
      %ws(j) = temp/n;
      
      %musum = zeros(1,d);
      %for i=1:n
      %  musum = musum + eta(i,j) * xs(i,:);
      %end
      %musum = musum + diff*0.01 * rand(1,d)-0.5*ones(1,d);
      %mus(j,:) = musum ./ temp;
      
    %  sigmasum = zeros(d,d);
    %  for i=1:n
    %    xmu = xs(i,:) - mus(j);
    %    sigmasum = sigmasum + eta(i,j) * xmu' * xmu;
    %  end
    %  sigmasum = sigmasum + diff*0.01 * rand(d,d)-0.5*ones(d,d);
    %  sigmas(j,:,:) = sigmasum ./ temp;
    end
    count = count + 1;
    %完了
    %eta
    %ws
    %mus
    %sigmas
    %L = - LogLH(xs,ws,mus,sigmas);
    %L
    %diff = L - lastL;
    %fprintf('count:%d  -> %10.5f\n',count,diff);
    %pause;
    %lastL = L;
  end
end


%データ列X,mu,sigmaを与えるとその正規分布における尤度を計算して返す
function y = LogLH(xs,ws,mus,sigmas)
  [n,d] = size(xs);
  prod = 1;
  for i = 1:n
    ax = xs(i,:)';
    prod = prod * q(ax,ws,mus,sigmas);
  end
  y = log(prod);
end

%モデルq
function y = q(ax,ws,mus,sigmas)
  [m d] = size(mus);
  y=0;
  for j=1:m
    y = y + ws(j) * N(ax,mus(j,:)',sigmas(j,:,:));
  end
end
%axは列ベクトル
function y = N(ax,mu,sigma)
  d = size(ax);
  d = d(1);
  y= 1 / ( (2 * pi)^(d/2) * det(sigma)^(1/2)) * exp(-1/2 * (ax-mu)' * pinv(sigma) * (ax-mu));
end