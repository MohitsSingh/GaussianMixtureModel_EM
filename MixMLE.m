
%�����K�E�X���f���ɑ΂��čŖސ�����s��
%�^����ꂽ�W�{�A�������ɑ΂��čœK�ȃp�����[�^��Ԃ�
%xs n*d�s�� �e�s��d�����W�{�������Ă���
%m ������
%e ��������Ɏg������ eta�̍��������菬�����Ȃ�����I������
%�Ԃ�l
%ws = m * 1 �̗�x�N�g��
%mus = m*d �̍s�� mu�̍s�x�N�g��������ł���
%sigmas = d*d*n cov��n������ł���
function [ws,mus,sigmas] = MixMLE(xs,m,e)
  [n,d] = size(xs);
  %�����l�ݒ�
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
    %E�X�e�b�v
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
    
    %�����v�Z
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
    
    %M�X�e�b�v
    %temp = m * 1 �s��
    temp = (sum(eta))';
    
    %�d�݂̌v�Z
    ws = temp ./ n;
    
    %���ς̌v�Z
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
    
    %y = n*m*d �s��
    temp5 = temp3 - temp4;
    
    temp6 = repmat(temp5,[1,1,1,d]);
    temp7 = permute(temp6,[1,2,4,3]);
    
    %Y = n*m*d*d�s��
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
    %����
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


%�f�[�^��X,mu,sigma��^����Ƃ��̐��K���z�ɂ�����ޓx���v�Z���ĕԂ�
function y = LogLH(xs,ws,mus,sigmas)
  [n,d] = size(xs);
  prod = 1;
  for i = 1:n
    ax = xs(i,:)';
    prod = prod * q(ax,ws,mus,sigmas);
  end
  y = log(prod);
end

%���f��q
function y = q(ax,ws,mus,sigmas)
  [m d] = size(mus);
  y=0;
  for j=1:m
    y = y + ws(j) * N(ax,mus(j,:)',sigmas(j,:,:));
  end
end
%ax�͗�x�N�g��
function y = N(ax,mu,sigma)
  d = size(ax);
  d = d(1);
  y= 1 / ( (2 * pi)^(d/2) * det(sigma)^(1/2)) * exp(-1/2 * (ax-mu)' * pinv(sigma) * (ax-mu));
end