%A-稀疏系数矩阵
%D-字典/测量矩阵（已知）
%X-测量值矩阵（已知）
%K-稀疏度
function A=OMP(D,F,X,L)
X=double(X);
[n,P]=size(X);
[n,K]=size(D);
%按列操作，分别求出每一列对应的最相关的基
for k=1:P
%a:每一列对应的相关基的系数
a=[];
%取二维信号的每一列信号
x=X(:,k);
%初始残差
residual=x;
%indx：索引集，L：测量次数
indx=zeros(L,1);
for j=1:L
%D转置与残差residual相乘，得到residual与每一列的内积值
residual=double(residual);
D=double(D);
proj=D'*residual;
%找内积值最大值的位置，即最相关基的position
pos=find(abs(proj)==max(abs(proj)));
%若最大值不止一个，取第一个
pos=pos(1);
%位置存入索引集的第j值
indx(j)=pos;
%indx(1:j)表示第一列前j个元素
%pinv：Pseudoinverse伪逆矩阵，一般用于处理长方形矩阵的求逆
%得到其相关基的对应系数，AD=X，A=inv(D)*X
%一般应该通过最小二乘来求
% a=pinv(D(:,indx(1:j)))*x;
An=D(:,indx(1:j));
a=inv(An' * An) * An'*x;
%继续求残差
residual=x-D(:,indx(1:j))*a;
end
%通过上面的循环，得到第k列的相关基对应的索引位置
temp=zeros(K,1);
temp(indx)=a;
%只显示非零值及其位置
%最终求得整个A，代入AD=X，即可求解
A(:,k)=temp;
end
% R=A'*D';
% R=uint8(R);
% imshow(R);