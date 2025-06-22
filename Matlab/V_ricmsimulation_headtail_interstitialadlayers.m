%% V: Head-Tail Simulation, Interstitial Adlayers
clear; clc;

%IRM Light Wavelength
lambda1 = 450e-9;
lambda2 = 532e-9; 
lambda3 = 610e-9;

%Refractive Indices
n0 = 1.520; %Substrate refractive index
n1 = 1.335; %Base adlayer refractive index
n2 = 1.490; %Phosphocholine lipid head group refractive index
n3 = 1.440; %Lipid tail group refractive index (octadecane used in model)
n4 = 1.335; %Interstitial adlayer refractive index
n5 = 1.335; %PBS top solution and water adlayer refractive index

%Distances and Wavelength
m = 0:30;
npoint=numel(m);
d1 = 1e-9; %Base adlayer thickness
d2 = 0.5e-9; %Phosphocholine head thickness
d3 = 4e-9; %Tail thickness
d4 = 0.5e-9; %Interstitial adlayer thickness

%Fresnel Reflectance Coefficients
r01 = (n0-n1)/(n0+n1); %Substrate to base adlayer
r12 = (n1-n2)/(n1+n2); %Base adlayer to head
r23 = (n2-n3)/(n2+n3); %Head to tail
r32 = -r23; %Tail to head
r24 = (n2-n4)/(n2+n4); %Head to interstitial adlayer
r42 = -r24; %Interstitial adlayer to head
r25 = (n2-n5)/(n2+n5); %Head to top solution

%Fresnel Transmittance Coefficients
t01 = 2*n0/(n0+n1);  %Substrate to base adlayer
t12 = 2*n1/(n1+n2); %Base adlayer to head
t23 = 2*n2/(n2+n3); %Head to tail
t32 = 2*n3/(n2+n3); %Tail to head
t24 = 2*n2/(n2+n4); %Head to interstitial adlayer
t42 = 2*n4/(n2+n4); %Interstitial adlayer to head
t25 = 2*n2/(n2+n5); %Head to top solution

%Transfer Matrices between layers
M01 = [1,r01;r01,1]/t01; %Substrate to base adlayer
M12 = [1,r12;r12,1]/t12; %Base adlayer to head
M23 = [1,r23;r23,1]/t23; %Head to tail
M32 = [1,r32;r32,1]/t32; %Tail to head
M24 = [1,r24;r24,1]/t24; %Head to interstitial adlayer
M42 = [1,r42;r42,1]/t42; %Interstitial adlayer to head
M25 = [1,r25;r25,1]/t25; %Head to top solution 

%Transfer Matrices within layers
phi1_l1 = 2*pi*n1*d1/lambda1; %Base adlayer
phi2_l1 = 2*pi*n2*d2/lambda1; %Head
phi3_l1 = 2*pi*n3*d3/lambda1; %Tail
phi4_l1 = 2*pi*n4*d4/lambda1; %Interstitial adlayer
M1_l1 = [exp(1i*phi1_l1),0;0,exp(-1*1i*phi1_l1)]; %Base adlayer
M2_l1 = [exp(1i*phi2_l1),0;0,exp(-1*1i*phi2_l1)]; %Head
M3_l1 = [exp(1i*phi3_l1),0;0,exp(-1*1i*phi3_l1)]; %Tail
M4_l1 = [exp(1i*phi4_l1),0;0,exp(-1*1i*phi4_l1)]; %Interstitial adlayer

phi1_l2 = 2*pi*n1*d1/lambda2; %Base adlayer
phi2_l2 = 2*pi*n2*d2/lambda2; %Head
phi3_l2 = 2*pi*n3*d3/lambda2; %Tail
phi4_l2 = 2*pi*n4*d4/lambda2; %Interstitial adlayer
M1_l2 = [exp(1i*phi1_l2),0;0,exp(-1*1i*phi1_l2)]; %Base adlayer
M2_l2 = [exp(1i*phi2_l2),0;0,exp(-1*1i*phi2_l2)]; %Head
M3_l2 = [exp(1i*phi3_l2),0;0,exp(-1*1i*phi3_l2)]; %Tail
M4_l2 = [exp(1i*phi4_l2),0;0,exp(-1*1i*phi4_l2)]; %Interstitial adlayer

phi1_l3 = 2*pi*n1*d1/lambda3; %Base adlayer
phi2_l3 = 2*pi*n2*d2/lambda3; %Head
phi3_l3 = 2*pi*n3*d3/lambda3; %Tail
phi4_l3 = 2*pi*n4*d4/lambda3; %Interstitial adlayer
M1_l3 = [exp(1i*phi1_l3),0;0,exp(-1*1i*phi1_l3)]; %Base adlayer
M2_l3 = [exp(1i*phi2_l3),0;0,exp(-1*1i*phi2_l3)]; %Head
M3_l3 = [exp(1i*phi3_l3),0;0,exp(-1*1i*phi3_l3)]; %Tail
M4_l3 = [exp(1i*phi4_l3),0;0,exp(-1*1i*phi4_l3)]; %Interstitial adlayer

%Light Pathway: m01 - m1 - m12 - (j-1)[m2 - m23 - m3 - m32 - m2 - m24 - m4
%- m42] - m2 - m23 - m3 - m32 - m2 - m25
for j=0:(npoint-1)
	M_l1 = M01*M1_l1*M12*[(M2_l1*M23*M3_l1*M32*M2_l1*M24*M4_l1*M42)^(j)]*M2_l1*M23*M3_l1*M32*M2_l1*M25
	r_l1(j+1)=M_l1(2,1)/M_l1(1,1);
    	t_l1(j+1)=1/M_l1(1,1);
	
	M_l2 = M01*M1_l2*M12*[(M2_l2*M23*M3_l2*M32*M2_l2*M24*M4_l2*M42)^(j)]*M2_l2*M23*M3_l2*M32*M2_l2*M25
	r_l2(j+1)=M_l2(2,1)/M_l2(1,1);
    	t_l2(j+1)=1/M_l2(1,1);
	
    	M_l3 = M01*M1_l3*M12*[(M2_l3*M23*M3_l3*M32*M2_l3*M24*M4_l3*M42)^(j)]*M2_l3*M23*M3_l3*M32*M2_l3*M25
	r_l3(j+1)=M_l3(2,1)/M_l3(1,1);
    	t_l3(j+1)=1/M_l3(1,1);
end
for k=1:npoint
        R_l1(k)=norm(r_l1(k))^2;
        T_l1(k)=norm(t_l1(k))^2*n5/n0;
		
		R_l2(k)=norm(r_l2(k))^2;
        T_l2(k)=norm(t_l2(k))^2*n5/n0;
		
		R_l3(k)=norm(r_l3(k))^2;
        T_l3(k)=norm(t_l3(k))^2*n5/n0;
end

Rbg=(((n0-n5)./(n0+n5)).^2);
CIRM_l1=R_l1./Rbg;
CIRM_l2=R_l2./Rbg;
CIRM_l3=R_l3./Rbg;


%Creating the graph
figure(2);
plot(m,CIRM_l1,'b--o')
hold on 
plot(m,CIRM_l2,'g--o')
plot(m,CIRM_l3,'r--o')
hold off
xlabel('Lipid Layers')
ylabel('IRM I/I_0')
grid on;
ylim([0,1.1])
xlim([min(m)-0.5,max(m)+0.5])
textct=sprintf('n_S_o_l=%4.3f\n n_I_A_d=%4.3f\n n_T_a_i_l%4.3f\n n_H_e_a_d=%4.3f\n n_B_A_d=%4.3f\n',n5,n4,n3,n2,n1,n0);
% text('Position',[20 1 0]);
% text(min(xlim())+0.6*(max(xlim())-min(xlim())),0.8*min(ylim())+0.2*(max(ylim())-min(ylim())),textct);
text(7,0.9,textct);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
opfnhead=sprintf('nSol=%4.2f,nIAd=%4.2f, nTail =%4.2f,nHead=%4.2f, nBAd = %4.2f, wavelength=%4.1fnm, HeadTailInterstitialAdlayers',n5,n4,n3,n2,n1,lambda1*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l1],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l1(2);
Mimg(412:612,412:612)=I_glass*CIRM_l1(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);

opfnhead=sprintf('nSol=%4.2f,nIAd=%4.2f, nTail =%4.2f,nHead=%4.2f, nBAd = %4.2f, wavelength=%4.1fnm, HeadTailInterstitialAdlayers',n5,n4,n3,n2,n1,lambda2*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l2],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l2(2);
Mimg(412:612,412:612)=I_glass*CIRM_l2(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);

opfnhead=sprintf('nSol=%4.2f,nIAd=%4.2f, nTail =%4.2f,nHead=%4.2f, nBAd = %4.2f, wavelength=%4.1fnm, HeadTailInterstitialAdlayers',n5,n4,n3,n2,n1,lambda3*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l3],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l3(2);
Mimg(412:612,412:612)=I_glass*CIRM_l3(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);
