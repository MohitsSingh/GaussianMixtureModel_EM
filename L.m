
%データ列X,mu,sigmaを与えるとその正規分布における尤度を計算して返す
function y = L(x,mu,sigma)
  [n,d] = size(x);
  invS = inv(sigma);
  prod = 1;
  %sum = sum + (n*d)/2 * log(2*pi) + n/2*log(det(sigma));
  for i = 1:n
    ax = (x(i,:) - mu)';
    %sum = sum + ax' * invS * ax;
    q = 1 / ( (2 * pi)^(d/2) * (det(sigma))^(1/2)) * exp(-1/2 * ax' * invS * ax) ;
    prod = prod * q;
    %sum = sum + log(q);
  end
  y = prod;
end