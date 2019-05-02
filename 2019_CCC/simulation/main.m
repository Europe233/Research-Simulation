clear;
close all;
%% Basic Matrixs
global n fre
N=6;% ½ÚµãÊý
L=[1 -1 0 0 0 0;
    -1 3 -1 -1 0 0;
    0 -1 2 0 -1 0;
    0 -1 0 2 -1 0;
    0 0 -1 -1 3 -1;
    0 0 0 0 -1 1];
[V,D]=eig(eye(N)+L);
n=3;
A=[ 1.8 0 -0.9;
    3.4 -0.2 -1.6;
    2.1 0.3 -1.1];
t=20;
fre=1000;
%alpha=2*norm(A)/D(1,1);
alpha=11;
gamma=20;
delta=0.01;
G1=diag([3,1,2,2,0,0]);
G2=[[1 0 1 1;
    0 1 0 0;
    1 0 0 1;
    0 1 1 0],zeros(4,2);
    zeros(2,6)];
X=zeros(n*N,t*fre+1);
X_dot=zeros(n*N);
Xi=zeros(n*N,t*fre+1);
Xi_dot=zeros(n*N);
%% Initial States
X(:,1)=(rand(n*N,1)-0.5)*8;
Xi(:,1)=(rand(n*N,1)-0.5)*10;
epsilon=zeros(n,t*fre+1);
time=linspace(0,t,fre*t+1);
%% Control
for i=1:t*fre
    H=H_matrix(i);
    Xi_dot=-gamma*(kron(L,eye(n))*X(:,i)+delta*Xi(:,i));
    X_dot=(kron(eye(N),A)-alpha*kron(L,eye(n)))*X(:,i)+kron(L,eye(n))*Xi(:,i)-alpha*kron(G1,eye(n))*X(:,i)+alpha*kron(G2,eye(n))*H;
    Xi(:,i+1)=Xi(:,i)+Xi_dot*1/fre;
    X(:,i+1)=X(:,i)+X_dot*1/fre;
end
%% Plot Results
for i=1:(t*fre+1)
    H=H_matrix(i);
    epsilon(:,i)=(kron(ones(1,N),eye(n))*kron(G2,eye(n))*H)/(ones(1,N)*G2*ones(N,1));
end
for j=1:n
    figure;
    for i=0:N-1
        plot(time,X(j+i*n,:),'LineWidth',1.5);
        hold on;
    end
    plot(time,epsilon(j,:),'k-.','LineWidth',2);
    hold on;
    xlabel('\bf{t} $\bf{(s)}$','FontSize',18,'Interpreter','latex');
    if(j==1)
        ylabel('$$x_i^{(1)}(t)$$','Interpreter','latex','FontSize',15);
        print -depsc -r600 x1
        handle=legend('$x_1^{(1)}(t)$','$x_2^{(1)}(t)$','$x_3^{(1)}(t)$','$x_4^{(1)}(t)$','$x_5^{(1)}(t)$','$x_6^{(1)}(t)$','$\epsilon^{(1)}(t)$');
    end
    if(j==2)
        ylabel('$$x_i^{(2)}(t)$$','Interpreter','latex','FontSize',15);
        print -depsc -r600 x2
        handle=legend('$x_1^{(2)}(t)$','$x_2^{(2)}(t)$','$x_3^{(2)}(t)$','$x_4^{(2)}(t)$','$x_5^{(2)}(t)$','$x_6^{(2)}(t)$','$\epsilon^{(2)}(t)$');
    end
    if(j==3)
        ylabel('$$x_i^{(3)}(t)$$','Interpreter','latex','FontSize',15);
        print -depsc -r600 x3
        handle=legend('$x_1^{(3)}(t)$','$x_2^{(3)}(t)$','$x_3^{(3)}(t)$','$x_4^{(3)}(t)$','$x_5^{(3)}(t)$','$x_6^{(3)}(t)$','$\epsilon^{(3)}(t)$');
    end
    
    set(handle, 'interpreter','latex');  
end