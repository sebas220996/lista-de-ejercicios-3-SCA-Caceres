// CONTROL SUSPENSION DE UN AUTOMOVIL
// Alumno: Sebastian Miguel Caceres Huaman
// Author: juan C. Cutipa-Luque
// Docente: juan C. Cutipa-Luque

//Codigo Basado en el  ejemplo del Prof. Rodriguez
//Adaptado  para el  modelo de suspension de un  automovil
clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows

///// PREGUNTA 1 ////////
// Parametros de la planta /////
k1 = 1;
k2 = 4;
b1 = 0.2;
b2 = 0.1;
m1 = 1;
m2 = 2;
// Modelo de la planta //
ap = [0 0 1 0; 0 0 0 1; -k1/m1 k1/m1 -b1/m1 b1/m1; k1/m2 -(k1 + k2)/m2 b2/m2 -(b1 + b2)/m2];
bp = [0 0; 0 0; 1/m1 0; 0 1/m2];
cp = [1 0 0 0; 0 1 0 0];
dp = 0*ones (2,2);
 
Gs= syslin("c",ap, bp, cp, dp)

//// PREGUNTA 2 //////
// Polos y Ceros de la planta //
[h1]=ss2tf(Gs)
p=spec(ap)
scf(1)
plzr(h1);
xtitle("Polos y zeros matrices de trnasmision de la planta")

////// PREGUNTA 3 /////
// Barreras de Estabilidad //
w=logspace(-3,3,400);
a=200;b=130;c=70;
x1=20*ones(1,a); x2=60*zeros(1,b); x3=-20*ones(1,c);
xt=[x1 x2 x3];
x4=[5*ones(1,100) 0*zeros(1,300)];

scf(2);
plot2d("ln", w, x4,3, rect=[10^-1 -60 10^3 60])
plot2d("ln", w, xt)
xgrid(12)
xtitle("Barreras de Estabilidad","Frequency w(rad/s)", "Amplitude (dB)");

////// PREGUNTA 4 //////
// Valores singulares de la planta //
w1 = logspace(-3,3);
sv = svplot(Gs,w1);
scf(3)
plot2d("ln", w1, [20*log(sv')/log(10)])
xgrid(12)
xtitle("Singular values plot planta G_s","Frequency (rad/s)", "Amplitude (dB)");

////// PREGUNTA 5 ///////
// Escalonamiento a la planta //
su = diag( [0.9614, 0.2753] )
sx = diag( [3.157, 11.47, 3.157, 11.47] )
sy = diag( [3.157 3.157] )
 
ap_ = sx*ap*inv(sx)
bp_ = sx*bp*inv(su)
cp_ = sy*cp*inv(sx)
dp_ = sy*dp*inv(su)

Gs_= syslin("c",ap_, bp_, cp_, dp_)

// Valores singulares de la planta escalonada //
sv1 = svplot(Gs_,w1);
scf(4)
plot2d("ln", w1, [20*log(sv1')/log(10)])
xgrid(12)
xtitle("Singular values plot planta escalonada G_s","Frequency (rad/s)", "Amplitude (dB)");

////// PREGUNTA 6 //////
// Controllability and Observability de la planta
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap,bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)

// Controllability and Observability de la planta escalonada
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc_ = cont_mat(ap_,bp_)
rankCc_=rank(Cc_)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O_ = obsv_mat(ap_, cp_)
rankO_=rank(O_)

// Planta aumentada con integradores antes del proyecto de controlador

[ns,nc] = size(bp_);                     //ns = número de entradas;  
                                        //nc = número de controles;   
a_1 = [ap_            bp_ ; 
       0*ones(nc,ns)  0*ones(nc,nc) ];
b_1 = [0*ones(ns,nc); eye(nc,nc)];
c_1 = [cp_      0*ones(nc,nc)];
d_1 = 0*ones(nc,nc)

Gs_1= syslin("c",a_1, b_1, c_1, d_1)

// Valores singulares de la planta escalonada con  el  integrador //
sv2 = svplot(Gs_1,w1);
scf(5)
plot2d("ln", w1, [20*log(sv2')/log(10)])
xgrid(12)
xtitle("Singular values plot planta escalonada con integradores G_s","Frequency (rad/s)", "Amplitude (dB)");

// LQR controller calculation
// Recuperar Target Loop resolviendo un problema de LQR barato
q = c_1'*c_1;          //Matriz de ponderación del estado
rho = 1e-9;        //Parámetro de recuperación de control barato

r = rho*eye(nc,nc)                 //Matriz de ponderación de control
         //how we calculate B
B=b_1*inv(r)*b_1';
A=a_1;
        //Solv the ricatti equation
X=riccati(A,B,q,'c','eigen');
        //matriz de ganacia
G_1=inv(r)*b_1'*X;

////// PREGUNTA 7 //////
//calculate observer Kalman Filter
ll =  inv(cp_*inv(-ap_)*bp_ + dp_);     
lh = -inv(ap_)*bp_*ll;
l = [lh                        //ll, lh - Para la conformación de 
     ll];                      //bucles de baja y alta frecuencia.

Gs_2= syslin("c",a_1, l, c_1, d_1)

// Valores singulares del filtro de bucle Abierto // 
sv3 = svplot(Gs_2,w1);
scf(6)
plot2d("ln", w1, [20*log(sv3')/log(10)])
xgrid(12)
xtitle("Singular values del filtro de bucle abierto","Frequency (rad/s)", "Amplitude (dB)");

// Filtro de Kalman
pnint=eye(nc,nc)               //Proceso  de matriz de intensidad de ruido
mu=0.01;                       //Medicion de la intensidad de ruido

mnint=mu*eye(nc,nc)            //Matriz de intensidad de ruido  de medicion

Ch=l*l';                       //Forma de Ch para "riccati" segun  Scilab
Ah=a_1';                       //Forma de Ah para "riccati" segun  Scilab

Bh=c_1'*inv(mnint)*c_1;

Xh=riccati(Ah,Bh,Ch,'c','eigen');

                           //ganacia H
H_1=(inv(mnint)*c_1*Xh)';

Gs_3= syslin("c",a_1, H_1, c_1, d_1)

// Valores singulares del observador Filtro  Kalman
sv4 = svplot(Gs_3,w1);
scf(7)
plot2d("ln", w1, [20*log(sv4')/log(10)])
xgrid(12)
xtitle("Singular values plot del filtro de kalman","Frequency (rad/s)", "Amplitude (dB)");


//ACTIVAR ESTA PARTE PARA LOS POLOS Y CEROS DEL FILTRO DE KALMAN GANACIA H 
//[h2]=ss2tf(Gs_3)
//scf(8)
//plzr(h2);
//xtitle("Polos y zeros del filtro ganancia H")

/////// PREGUNTA 8 ////////
// COMPENSADOR K(S) DE LA FORMA DEL PORF. RODRIGUEZ //
ak = [ a_1-b_1*G_1-H_1*c_1  0*ones(ns+nc,nc)
       G_1                  0*ones(nc,nc) ];
bk = [ H_1
       0*ones(nc,nc) ];
ck = [0*ones(nc, ns+nc) eye(nc,nc) ];
dk = [0*eye(nc,nc)];

Gs_4= syslin("c",ak, bk, ck, dk)
// Valores singulares del compensador "K(s)" //
sv4 = svplot(Gs_4,w1);
scf(9)
plot2d("ln", w1, [20*log(sv4')/log(10)])
xgrid(12)
xtitle("Singular values plot del compensador Ks","Frequency (rad/s)", "Amplitude (dB)");

/////// PREGUNTA 9 ///////
//SENSIBILIDAD "S" Y SENSIBILIDAD COMPLEMENTARIA "T"
//Analisis en bucle abierto
al = [ ap_                     bp_*ck
       0*ones(ns+nc+nc,ns)     ak    ];
bl = [ 0*ones(ns,nc)
       bk ];
cl = [ cp_  0*ones(nc,ns+nc+nc) ];
dl = [0*eye(nc,nc)];

Gs_5= syslin("c",al, bl, cl, dl)

// Valores Singulares de bucle abierto//
sv5 = svplot(Gs_5,w1);
scf(10)
plot2d("ln", w1, [20*log(sv5')/log(10)])
xgrid(12)
xtitle("Singular values plot del bucle abierto","Frequency (rad/s)", "Amplitude (dB)");

// Valores singulares de sensibilidad S //
Gs_6= syslin("c",al-bl*cl, bl, -cl, eye(nc,nc))
sv6 = svplot(Gs_6,w1);
scf(11)
plot2d("ln", w1, [20*log(sv6')/log(10)])
xgrid(12)
xtitle("Singular values plot Sensibilidad","Frequency (rad/s)", "Amplitude (dB)");

// Valores singulares de sensibilidad complementaria T //
Gs_7= syslin('c',al-bl*cl, bl, cl, dl)
sv7 = svplot(Gs_7,w1);
scf(12)
plot2d("ln", w1, [20*log(sv7')/log(10)])
xgrid(12)
xtitle("Singular values plot Sensibilidad Complementaria","Frequency (rad/s)", "Amplitude (dB)");

// Valores Singulares de S y T juntos // 
scf(13)
plot2d("ln", w1, [20*log(sv6')/log(10) 20*log(sv7')/log(10)],leg="Sensibilidad@Sensibilidad Complementaria");
xgrid(12)
xtitle("Singular values plot Sensibilidad y Sensi. Complementaria","Frequency (rad/s)", "Amplitude (dB)");

//RESPUESTA AL STEP DEL SISTEMA CONTROLADO
//la respuesta al step se realizo en  el  programa GNU OCTAVE ONLINE pasando los datos 
//necesarios para la grafica de esta misma
