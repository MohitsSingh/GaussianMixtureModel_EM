
%d*n*mを返す
function xx = generatedata(d,n)
  
  [m,qm] = settings();
  
  %データはd次元の列ベクトル
  %全てのデータ入れ
  xx = zeros(d,n,m);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %   訓練データ生成
  %%%%%%%%%%%%%%%%%%%%%%%%%

  %1共分散!=0
  A_3 = randn(d,d);
  A_3 = A_3' * A_3;
  mu = randn(d,1) * ones(1,n);
  xx(:,:,1) = sqrtm(A_3)*randn(d,n) + mu;
  
  %2共分散0
  A_2 = zeros(d,d);
  for i = 1:d
     A_2(i,i) = rand;
  end
  mu = randn(d,1) * ones(1,n);
  xx(:,:,2) = sqrtm(A_2)*randn(d,n)+ mu;

  %3共分散0分散同一
  A_1 = zeros(d,d);
  A_1(1,1) = rand;
  for i = 2:d
     A_1(i,i) = A_1(1,1);
  end
  mu = randn(d,1) * ones(1,n);
  xx(:,:,3) = sqrtm(A_1)*randn(d,n)+ mu;

  %4
  %[0,1]の一様分布
  xx(:,:,4) = rand(d,n);
  
  %5
  %[0,1]と[-1,0]の一様分布
  for count=1:n
    xx(:,count,5) = (randi(2)*2-3) * rand(d,1);
  end
  
  return;
  
  %6
  %0.5を中心とした円
  count = 1;
  while 1
    if count == n
      break
    end
    u = rand(d,1);
    x = (u - 0.5 * ones(d,1)) * 2;
    r = x' * x;
    if 1 < r 
      continue;
    end
    xx(:,count,6) = u;
    count = count + 1;
  end
  
  %7
  %[0,k]の一様分布
  for count=1:n
    for i = 1:d
      xx(i,count,7) = randi(i+1) * rand(1,1);
    end
  end

end
