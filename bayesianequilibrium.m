function [A,payoff,iterations,err] = npg(M,U)
M = [2 2 2 2];
pr=[2/5 3/10 1/5 1/10];
K=[13/6 5/2 2 2 13/6 5/2 2 2 4/3 3/2 2 2 4/3 3/2 2 2 ;23/12 1/2 9/4 1/2 2/3 1/2 5/4 1/2 23/12  1/2 9/4 1/2 2/3 1/2 5/4 1/2];
U(1,1)=pr(1)*K(1,1)+pr(2)*K(1,5);
U(1,2)=pr(3)*K(1,9)+pr(4)*K(1,13);
U(1,3)=pr(1)*K(2,1)+pr(3)*K(2,9);
U(1,4)=pr(2)*K(2,5)+pr(4)*K(2,13);

U(2,1)=pr(1)*K(1,1)+pr(2)*K(1,6);
U(2,2)=pr(3)*K(1,9)+pr(4)*K(1,14);
U(2,3)=U(1,3);
U(2,4)=pr(2)*K(2,6)+pr(4)*K(2,14);

U(3,1)=pr(1)*K(1,2)+pr(2)*K(1,5);
U(3,2)=pr(3)*K(1,10)+pr(4)*K(1,13);
U(3,3)=pr(1)*K(2,2)+pr(3)*K(2,10);
U(3,4)=U(1,4);

U(4,1)=pr(1)*K(1,2)+pr(2)*K(1,6);
U(4,2)=pr(3)*K(1,10)+pr(4)*K(1,14);
U(4,3)=U(3,3);
U(4,4)=U(2,4);

U(5,2)=pr(3)*K(1,11)+pr(4)*K(1,15);
U(5,3)=pr(1)*K(2,1)+pr(3)*K(2,11);
U(5,4)=pr(2)*K(2,5)+pr(4)*K(2,15);

U(6,2)=pr(3)*K(1,11)+pr(4)*K(1,16);
U(6,3)=U(5,3);
U(6,4)=pr(2)*K(2,5)+pr(4)*K(2,16);

for i=5:8
    U(i,1)=U(i-4,1);
end

for i=7:8
    U(i,4)=U(i-2,4);
end

U(7,2)=pr(3)*K(1,12)+pr(4)*K(1,15);
U(7,3)=pr(1)*K(2,2)+pr(3)*K(2,12);
U(8,3)=U(7,3);
U(8,2)=pr(3)*K(1,12)+pr(4)*K(1,16);

for i=9:10
    U(i,1)=pr(1)*K(1,3)+pr(2)*K(1,i-2);
end
for i=11:12
    U(i,1)=pr(1)*K(1,4)+pr(2)*K(1,i-4);
end

for i=13:16
    U(i,1)=U(i-4,1);
end

for i=9:16
    U(i,2)=U(i-8,2);
end

U(9,3)=pr(1)*K(2,3)+pr(3)*K(2,9);
U(10,3)=U(9,3);
U(11,3)=pr(1)*K(2,4)+pr(3)*K(2,10);
U(12,3)=U(11,3);
U(13,3)=pr(1)*K(2,3)+pr(3)*K(2,11);
U(14,3)=U(13,3);
U(15,3)=pr(1)*K(2,4)+pr(3)*K(2,12);
U(16,3)=U(15,3);

U(9,4)=pr(2)*K(2,7)+pr(4)*K(2,13);
U(11,4)=U(9,4);
U(10,4)=pr(2)*K(2,8)+pr(4)*K(2,14);
U(12,4)=U(10,4);

U(13,4)=pr(2)*K(2,7)+pr(4)*K(2,15);
U(15,4)=U(13,4);
U(14,4)=pr(2)*K(2,8)+pr(4)*K(2,16);
U(16,4)=U(14,4);

p = 1; V = 1;
n = length(M);
s = sum(M);

A = zeros(max(M),n);
payoff = zeros(1,n);
for i = 1 : n
    p = p * M(1,i);
end
if p ~= size(U,1) || n ~= size(U,2)
    error('Error: Dimension mismatch!');
end
P = zeros(1,n);
N = zeros(1,n);
P(n) = 1;
for i = n-1 : -1 : 1
    P(i) = P(i+1) * M(i+1);
end
for i = 1 : n
    N(i) = p / P(i);
end
x0 = zeros(s,1); k = 1;
for i = 1 : n
    for j = 1 : M(1,i)
        x0(k) = 1 / M(1,i); k = k + 1;
    end
end
Us = sum(U,2);
for i = 1 : n
    V = V * (1 / M(i)) ^ M(i);
end
x0 = [x0 ; V * (sum(U,1)')];
Aeq = zeros(n,s+n); cnt = 0;
for i = 1 : n
    if i ~= 1
        cnt = cnt + M(i-1);
    end
    for j = 1 : s
        if j <= sum(M(1:i)) &&  j > cnt
            Aeq(i,j) = 1;
        end
    end
end
beq = ones(n,1);
I = ones(p,n);
counter = 0; count = 1;
for i = 1 : n
    for j = 1 : N(i)
        counter = counter + 1;
        if i ~= 1
            if counter > sum(M(1:i))
                counter = counter-M(i);
            end
        end
        for k = 1 : P(i)
            I(count) = counter;
            count = count + 1;
        end
    end
end
lb = zeros(s+n,1);
ub = ones(s+n,1);
pay = zeros(s,1);
counter = 0;
for i = 1 : n
    for j = 1 : M(i)
        counter = counter + 1;
        pay(counter) = i + s;
    end
end
for i = 1 : n
    lb(s+i) = -inf;
    ub(s+i) = inf;
end
[x,fval,exitflag,output] = gamer(n,Us,p,I,s,ub,lb,x0,Aeq,beq,pay,U);
count = 0;
for i = 1 : n
    for j = 1 : M(i)
        count = count + 1;
        A(j,i) = abs(x(count));
    end
    payoff(1,i) = x(s+i);
end
iterations = output.iterations;
err = abs(fval);