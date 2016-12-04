
%データ列X,mu,sigmaを与えるとその正規分布における尤度を計算して返す
function y = LogL(x,mu,sigma)
  [n,d] = size(x);
  invS = inv(sigma);
  
  sum = 0;
  sum = sum + (n*d)/2 * log(2*pi) + n/2*log(det(sigma));
  for i = 1:n
    ax = (x(i,:) - mu)';
    sum = sum + ax' * invS * ax;
  end
  y = - sum;
end