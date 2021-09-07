function code_ParallelServers_AoI_Energy


% uncomment the next line to generate Fig.3  
%do_plot_Fig3_Fig4_Fig5(1000)

% uncomment the next line to generate Fig.4 
%do_plot_Fig3_Fig4_Fig5(1)

% uncomment the next line to generate Fig.5 
%do_plot_Fig3_Fig4_Fig5(0.1)

% uncomment the next line to generate Fig.6 
%do_plot_Fig6_Fig7(0.1)

% uncomment the next line to generate Fig.7 
%do_plot_Fig6_Fig7(0.45)

end


function do_plot_Fig6_Fig7(p)

mu=1;
v_a=[.1:.1:1 1:100 100:20:1000];
ind=1;
%p=0.1;
for a=v_a
    l=10; r1(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    l=20; r2(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    l=50; r3(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    l=100; r4(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    l=200; r5(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    l=500; r6(ind)=compute_aaoi_twoparallel(l*p,l*(1-p),a,mu);
    ind=ind+1;
end
figure; 
hold on
semilogx(v_a,r1)
semilogx(v_a,r2)
semilogx(v_a,r3)
semilogx(v_a,r4)
semilogx(v_a,r5)
semilogx(v_a,r6)
legend('\lambda=10','\lambda=20','\lambda=50', '\lambda=100',...
    '\lambda=200','\lambda=500')

end


function do_plot_Fig3_Fig4_Fig5(a)

mu=1;
v_l=[.1:.1:1 1:100 100:20:1000];
ind=1;
for l=v_l
    r1(ind)=compute_aaoi_twoparallel(l/2,l/2,a,mu);
    r2(ind)=compute_aaoi_single(l/2,a,mu);
    r3(ind)=compute_aaoi_single(l,a,2*mu);
    ind=ind+1;
end
figure; 
hold on
semilogx(v_l,r1)
semilogx(v_l,r2)
semilogx(v_l,r3)
r1
r2
r3
legend('2 parallel', 'single l/2,mu','single l,2mu')

r1(end-30:end)
r2(end-30:end)
r3(end-30:end)
end




function aaoi=compute_aaoi_twoparallel(l1,l2,a,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stationnary distribution 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=[1, 1, 1, 1, 1, 1, 1, 1;... 
 -l1 - l2 - a, 0, mu, mu, mu, mu, 0, 0;...
  a, -l1 - l2, 0, 0, 0, 0, 0, 0;
  l2, 0, -(l1 + mu + a), 0, 0, 0, mu, mu;
  0, l2, a, -(l1 + mu), 0, 0, 0, 0;
  l1, 0, 0, 0, -(l2 + a + mu), 0, mu, mu; 
  0, l1, 0, 0, a, -(l2 + mu), 0, 0;
  0, 0, l1, 0, l2, 0, -(a + mu + mu), 0];
n=[1,0,0,0,0,0,0,0]';
vpi=M\n;
pi000=vpi(1);pi001=vpi(2);pi010=vpi(3);pi011=vpi(4);
pi100=vpi(5);pi101=vpi(6);pi110=vpi(7);pi111=vpi(8);
l=l1+l2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Average AoI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 v1=[-(l+a),0,mu,0,0,mu,mu,0,0,mu,0,0,0,0,0,0];
 v2=[a,-l,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
 v3=[l2, 0, -(l1+a+mu),0,0,0,0,0,0,0,0,0,0,0,mu,0];
 v4=[0, 0,0,-(l+a+mu),0,0,0,0,0,0,0,0,0,0,0,mu];
 v5=[0,l2,a,0,-(l1+mu),0,0,0,0,0,0,0,0,0,0,0];
 v6=[0,0,0,a,0,-(l+mu),0,0,0,0,0,0,0,0,0,0];
 v7=[l1,0,0,0,0,0,-(l2+a+mu),0,0,0,0,0,0,0,0,mu];
 v8=[0,0,0,0,0,0,0,-(l+a+mu),0,0,0,0,0,0,mu,0];
 v9=[0,l1,0,0,0,0,0,0,-(l2+mu),0,0,0,0,0,0,0];
 v10=[0,0,0,0,0,0,0,a,0,-(l2+mu),0,0,0,0,0,0];
 v11=[0,0,l1,0,0,0,l2,0,0,0,-(2*mu+a),0,0,0,0,0];
 v12=[0,0,0,0,0,0,0,l2,0,0,0,-(l1+2*mu+a),0,0,0,0];
 v13=[0,0,0,l1,0,0,0,0,0,0,0,0,-(l2+2*mu+a),0,0,0];
 v14=[0,0,0,0,l1,0,0,0,l2,0,a,0,0,-2*mu,0,0];
 v15=[0,0,0,0,0,0,0,0,0,l2,0,a,0,0,-(l1+2*mu),0];
 v16=[0,0,0,0,0,l1,0,0,0,0,0,0,a,0,0,-(l2+2*mu)];
 A=-[v1;v2;v3;v4;v5;v6;v7;v8;v9;v10;v11;v12;v13;v14;v15;v16];
b=[vpi(1) vpi(2) vpi(3) vpi(3) vpi(4) vpi(4) vpi(5) vpi(5) vpi(6) vpi(6) vpi(7) vpi(7) vpi(7) vpi(8) vpi(8) vpi(8)]';
x=A\b;
aaoi=x(1)+x(2)+x(3)+x(5)+x(7)+x(9)+x(11)+x(14);

      
end



function aaoi=compute_aaoi_single(l,a,mu)
b=0;
vpi=compute_stat_dist_single(a,b,l,mu);
A=(vpi(4)+(a*vpi(2))/(l+mu+a))/(l+mu+b-(a*b)/(l+mu+a));
N=[l+a, -mu,-b, 0;
   -l ,mu+a, 0,-b;
   -a ,  0,l+b, 0;
   0,   -a, -l,mu+b];
v=[vpi(1)+mu*A,vpi(2),vpi(3),vpi(4)]';
aaoi=sum(N\v);
end

function pi_=compute_stat_dist_single(a,b,l,mu)

b=0;
M=[a+l,-mu , -b,-mu; 
   -l ,mu+a,  0,-b; 
   -a ,   0,l+b, 0;
   1  ,   1,  1, 1];
n=[0,0,0,1]'; pi_=M\n;


end
