%�K�E�X�������f���ɑ΂��čŖސ�����s���œK�ȃp�����[�^��Ԃ��֐�
%����
% xs d*n�s�� �e�s��d�����W�{�������Ă���
% m ������
% e ��������Ɏg������ eta�̍��������菬�����Ȃ�����I������
%�Ԃ�l
% ws = m * 1 �̗�x�N�g��
% mus = d*m �̍s�� mu�̍s�x�N�g��������ł���
% sigmas = d*d*n ���U�����U�s��n������ł���
function [ws,mus,sigmas] = MixMLE1(xs,m,e)
  [d,n] = size(xs);
  
  %�����l�ݒ�
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
  
  while 1
    lasteta = eta;
    %E�X�e�b�v
    for i=1:n
      x = xs(:,i);
      temp = 0;
      for j=1:m
        eta(i,j) = ws(j) * N(x,mus(:,j),sigmas(:,:,j));
        temp = temp + eta(i,j);
      end
      eta(i,:) = eta(i,:) ./ temp;
    end
    
    %��������
    diff = sum(sum(~((eta - lasteta) < e)))
    if diff == 0
      break;
    end
    
    %M�X�e�b�v
    temp = (sum(eta))';
    %�d�݂̌v�Z
    ws = temp ./ n;
    %���ς̌v�Z
    temp1 = repmat(eta,[1,1,d]);%n,m,d
    temp1 = permute(temp1,[3,2,1]);%
    temp2 = repmat(xs,[1,1,m]);%d,n,m
    temp2 = permute(temp2,[1,3,2]);%d,m,n
    mus = sum((temp1.*temp2),3);%d,m
    mus = mus ./ (repmat(temp,[1,d]))';
    %���U�����U�s��̌v�Z
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
end

%���f��q�ɂ�����m�����v�Z���ĕԂ��֐�
%ax = d*1
%mu = d*m
%sigma = d*d*m
function y = q(ax,ws,mus,sigmas)
  [m d] = size(mus);
  y=0;
  for j=1:m
    y = y + ws(j) * N(ax,mus(:,j),sigmas(:,:,j));
  end
end

%���K���z�̊m����Ԃ��֐�
%ax = d*1
%mu = d*1
%sigma = d*d
function y = N(ax,mu,sigma)
  d = size(ax);
  d = d(1);
  y= 1 / ( (2 * pi)^(d/2) * det(sigma)^(1/2)) * exp(-1/2 * (ax-mu)' * pinv(sigma) * (ax-mu));
end
