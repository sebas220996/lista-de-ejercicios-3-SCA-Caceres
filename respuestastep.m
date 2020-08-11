clear;
clc;

k1 = 1;
k2 = 4;
b1 = 0.2;
b2 = 0.1;
m1 = 1;
m2 = 2;

ap = [0 0 1 0; 0 0 0 1; -k1/m1 k1/m1 -b1/m1 b1/m1; k1/m2 -(k1 + k2)/m2 b2/m2 -(b1 + b2)/m2];
bp = [0 0; 0 0; 1/m1 0; 0 1/m2];
cp = [1 0 0 0; 0 1 0 0];
dp = 0*ones (2,2);

su = diag( [0.9614, 0.2753] );
sx = diag( [3.157, 11.47, 3.157, 11.47] );
sy = diag( [3.157 3.157] );

ap_ = sx*ap*inv(sx);
bp_ = sx*bp*inv(su);
cp_ = sy*cp*inv(sx);
dp_ = sy*dp*inv(su);


[ns,nc] = size(bp_);

G_1=[31594.16384683239  -9.697152391243421   1339.452219279039   0.992351085342597  93.79144595181171   0.157413956229137;
    73.53305501875495   8690.277434064428   3.618587662394658   306.5880020165346  0.157413956229048   113.0201188405263];

H_1=[
   10.                 0.000000000000002
  -0.000000000000019   36.33196072220461
  -0.000000000000049   0.000000000000011
  -0.000000000000027  -0.000000000000017
   3.045296167247377  -3.045296167247384
  -0.872030408615777   4.360152043078865];


a_1 = [ap_            bp_ ; 
       0*ones(nc,ns)  0*ones(nc,nc) ];
b_1 = [0*ones(ns,nc); eye(nc)];
c_1 = [cp_      0*ones(nc,nc)];

ak = [ a_1-b_1*G_1-H_1*c_1  0*ones(ns+nc,nc)
       G_1                  0*ones(nc,nc) ];
bk = [ H_1
       0*ones(nc,nc) ];
ck = [0*ones(nc, ns+nc) eye(nc) ];

al = [ ap_                     bp_*ck
       0*ones(ns+nc+nc,ns)     ak    ];
bl = [ 0*ones(ns,nc)
       bk ];
cl = [ cp_  0*ones(nc,ns+nc+nc) ];
dl = [0*eye(nc)];


[y,t] = step(ss(al-bl*cl, bl, cl, 0*eye(nc)));
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
grid on
pause

plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
grid on
pause

plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
grid on
pause

plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
grid on
pause