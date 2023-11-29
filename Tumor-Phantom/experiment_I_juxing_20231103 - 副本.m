
clear
tic     

% ָ��Ҫ��ȡ�� Excel �ļ����ļ���
filename = '16chan Circular Scan Data-20231102-203027.xls';  % �滻Ϊʵ���ļ���

% ʹ�� xlsread ��������ȡ Excel �ļ��е�����
% ��һ��������������ݣ��ڶ�������������ı����ݣ����������������ԭʼ����
[data, text, raw] = xlsread(filename);

% ������ѡ�����ȡ���ݣ�������ȡ�ı����ݺ�ԭʼ����
% data = xlsread(filename);

F=60e6;

[cc,nn]=size(data);

T=1/F;
t_data=(1:cc)*T;


%%
%�ɼ�����
 for i=1:nn                 %��ȡȥ��ֱ������������
    
    S(:,i)=data(:,i);
    L = length(S(:,i));             % Number of signal ����һ������2����ָ����
                          %�������Զ�������ֵ�����Ի��С����������ʹ��ʵ�ʵķ����׼�С
Y = fft(S(:,i),L);   %����N����2������ʵ��ֵ��NԽ�󣬷�ֵ����Խ��
Y(1,1)=0;
V1=ifft(Y);
%V1((1350:end),1)=0;
V(:,i)=V1;
dVV(:,i)=diff(V(:,i))*F;
inputsignal=V(:,i);
cutoff_freq=[0.3e6 3e6];
[filtered_signal,filtb,filta]=bandpass_butterworth(inputsignal,cutoff_freq,F,3);
VV(:,i)=filtered_signal;
 end
%%
t_step=T;                            %���ݲɼ�ʱ�䲽��
cs=1500;                             %�趨��ѹ�����ٶ�
%offset=1680;
rd=0.1;                         %̽ͷ������֮��ľ���(m)
l_max=0.055;                         %���ͼ�񱳾��ĳ���m��
l_max1=0.055;                        %���ͼ�񱳾��Ŀ�m��

p_delta=VV;

%%  ɨ�����ǰ��׼��
[m,n]=size(p_delta);                %����ԭʼ���ݾ����С
Num_detector=n;                     %̽ͷ����
NumberOfPointschang=500;            %���ͼ�񱳾����ķֱ���
NumberOfPointskuan=500;             %���ͼ�񱳾���ķֱ���
[X,Y] = meshgrid(-l_max/2:l_max/(NumberOfPointschang-1):l_max/2,-l_max1/2:l_max1/(NumberOfPointskuan-1):l_max1/2);  %��Ŀ��������Ϊԭ�㽨����άֱ������ϵ
pt2=p_delta;                        %��ԭʼ���ݴ��ݳ�������ֹ�ƻ�ԭʼ����
%pt3=diff(pt2)/t_step;               %����ٶ��ƺ���
Z=zeros(NumberOfPointschang,NumberOfPointskuan);                   %Ԥ�������ڴ�λ��

%%  ��ʼ����ɨ��
for indexi=1:NumberOfPointschang                                  %��������ɨ��
    for indexj=1:NumberOfPointskuan                               %���������ɨ��
      pd=pt2;
      x0=X(1,indexi);
      y0=Y(indexj,1);
    t_step=1/F;                     %���ݲɼ�ʱ�䲽��
    [~,n_pt2]=size(pd);                             %����ɨ�����
    JBt=0;                                          %��ʼ�����������׼�����м���
    for i=1:n_pt2                                   %��һɨ������̽ͷ����
        alpha=-18.334/16;                            %����̽ͷÿ����ת�Ƕȣ��ȣ�
         rot_angle=-20/180*pi;
       t_delay=3.88e-5;    %3.32e-5 3.95e-5
        alpha_rd=alpha*(i-1)/180*pi+rot_angle;      %��̽ͷ�Ƕ�ת��Ϊ�����ƣ����ڼ���                
        rd_x=rd*cos(alpha_rd);                      %������ɨ��̽ͷ������ת��Ϊֱ�����꣨X������
        rd_y=rd*sin(alpha_rd);                      %������ɨ��̽ͷ������ת��Ϊֱ�����꣨Y������
        dst=(rd_x-x0)^2+(rd_y-y0)^2;                %ȷ��̽ͷ��ƽ��ֱ������ϵ��λ��
       tm=fix(sqrt(dst)/cs/t_step)-ceil(t_delay/t_step);                %�ҵ���·�ź��ڸ����ص��Ӧ����ѹֵ        
     %   d_nor=(rd-(x0*cos(alpha)+y0*sin(alpha)));   %��ⳤ����Դ��λ��
    %    JBt=JBt+d_nor*pd(tm,i)/dst;                 %������Դ�ۼ�
     if tm>m
      JBt=JBt;
  else
    JBt=JBt+pd(tm,i)/sqrt(dst);  
  end
            
    end
  %  JB=JBt*180/n_pt2*pi/180/(2*pi*cs^3);            %���������շֲ�         
  JB=JBt/(2*pi*cs^3);  
      
      Z(indexi,indexj)=JB;
    end
end
%% �������ͼ��
figure

set(gcf,'Position',[0 0 400 400]);
axis normal
% set(gca,'XLim',[0 12]);
% set(gca,'YLim',[0 12]);
% set(gca,'LooseInset',get(gca,'TightInset'));
% set(gca,'LooseInset',[0.5 0.5 0.5 0.5]);
image_nn=Z;
imagesc(image_nn);
colormap(gray)

