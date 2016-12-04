
%混合ガウスモデルに対して最尤推定を行う
%与えられた標本、混合数に対して最適なパラメータを返す
%xs n行列 各行にd次元標本が入っている
%m 混合数
%e 収束判定に使う正数 etaの差がこれより小さくなったら終了する
function [ws,mus,sigmas] = MixMLE1(xs,m,e)
  %dは1次元限定
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
  mus = rand(m);
  sigmas = rand(m);
  %pause;
  
  eta = zeros(n,m);
  lasteta = zeros(n,m);
  
  count = 0;
  lastL = 0;
  while 1
    lasteta = eta;
    %Eステップ
    for i=1:n
      x = xs(i);
      temp = 0;
      for j=1:m
        eta(i,j) = ws(j) * normpdf(x,mus(j),sigmas(j));
        temp = temp + eta(i,j);
      end
      eta(i,:) = eta(i,:) ./ temp;
    end
    
    %差分計算
    diff = sum(sum(~((eta - lasteta) < 0.001)))
    if diff == 0
      return
    end
    
    %Mステップ
    temp = sum(eta);
    ws = temp'./n;
    mu = (x'*eta ./ temp)';
    sigma = (x'*eta ./ temp)';
    
    for j=1:m
      temp = sum(eta(:,j));
      ws(j) = temp/n;
      
      musum = zeros(1,d);
      for i=1:n
        musum = musum + eta(i,j) * xs(i,:);
      end
      musum = musum + diff*0.01 * rand(1,d)-0.5*ones(1,d);
      mus(j,:) = musum ./ temp;
      repmat
      sigmasum = zeros(d,d);
      for i=1:n
        xmu = xs(i,:) - mus(j);
        sigmasum = sigmasum + eta(i,j) * xmu' * xmu;
      end
      sigmasum = sigmasum + diff*0.01 * rand(d,d)-0.5*ones(d,d);
      sigmas(j,:,:) = sigmasum ./ temp;
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