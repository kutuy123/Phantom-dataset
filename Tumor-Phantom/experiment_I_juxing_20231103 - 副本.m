
clear
tic     

% 指定要读取的 Excel 文件的文件名
filename = '16chan Circular Scan Data-20231102-203027.xls';  % 替换为实际文件名

% 使用 xlsread 函数来读取 Excel 文件中的数据
% 第一个输出参数是数据，第二个输出参数是文本数据，第三个输出参数是原始数据
[data, text, raw] = xlsread(filename);

% 您可以选择仅读取数据，而不读取文本数据和原始数据
% data = xlsread(filename);

F=60e6;

[cc,nn]=size(data);

T=1/F;
t_data=(1:cc)*T;


%%
%采集数据
 for i=1:nn                 %获取去除直流分量的数据
    
    S(:,i)=data(:,i);
    L = length(S(:,i));             % Number of signal 长度一般设置2的幂指数，
                          %不够会自动补充零值，所以会减小整体能量，使得实际的幅度谱减小
Y = fft(S(:,i),L);   %除以N乘以2才是真实幅值，N越大，幅值精度越高
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
t_step=T;                            %数据采集时间步长
cs=1500;                             %设定声压传播速度
%offset=1680;
rd=0.1;                         %探头与样本之间的距离(m)
l_max=0.055;                         %输出图像背景的长（m）
l_max1=0.055;                        %输出图像背景的宽（m）

p_delta=VV;

%%  扫描成像前的准备
[m,n]=size(p_delta);                %测量原始数据矩阵大小
Num_detector=n;                     %探头个数
NumberOfPointschang=500;            %输出图像背景长的分辨率
NumberOfPointskuan=500;             %输出图像背景宽的分辨率
[X,Y] = meshgrid(-l_max/2:l_max/(NumberOfPointschang-1):l_max/2,-l_max1/2:l_max1/(NumberOfPointskuan-1):l_max1/2);  %以目标体中心为原点建立二维直角坐标系
pt2=p_delta;                        %将原始数据传递出来，防止破坏原始数据
%pt3=diff(pt2)/t_step;               %求解速度势函数
Z=zeros(NumberOfPointschang,NumberOfPointskuan);                   %预留成像内存位置

%%  开始成像扫描
for indexi=1:NumberOfPointschang                                  %成像区域长扫描
    for indexj=1:NumberOfPointskuan                               %成像区域宽扫描
      pd=pt2;
      x0=X(1,indexi);
      y0=Y(indexj,1);
    t_step=1/F;                     %数据采集时间步长
    [~,n_pt2]=size(pd);                             %计算扫描次数
    JBt=0;                                          %初始化计算变量，准备进行计算
    for i=1:n_pt2                                   %逐一扫描所有探头数据
        alpha=-18.334/16;                            %计算探头每次旋转角度（度）
         rot_angle=-20/180*pi;
       t_delay=3.88e-5;    %3.32e-5 3.95e-5
        alpha_rd=alpha*(i-1)/180*pi+rot_angle;      %将探头角度转化为弧度制，便于计算                
        rd_x=rd*cos(alpha_rd);                      %将环形扫描探头极坐标转换为直角坐标（X分量）
        rd_y=rd*sin(alpha_rd);                      %将环形扫描探头极坐标转换为直角坐标（Y分量）
        dst=(rd_x-x0)^2+(rd_y-y0)^2;                %确定探头在平面直角坐标系的位置
       tm=fix(sqrt(dst)/cs/t_step)-ceil(t_delay/t_step);                %找到该路信号在该像素点对应的声压值        
     %   d_nor=(rd-(x0*cos(alpha)+y0*sin(alpha)));   %求解长点与源点位置
    %    JBt=JBt+d_nor*pd(tm,i)/dst;                 %进行声源累加
     if tm>m
      JBt=JBt;
  else
    JBt=JBt+pd(tm,i)/sqrt(dst);  
  end
            
    end
  %  JB=JBt*180/n_pt2*pi/180/(2*pi*cs^3);            %计算热吸收分布         
  JB=JBt/(2*pi*cs^3);  
      
      Z(indexi,indexj)=JB;
    end
end
%% 输出反演图像
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

