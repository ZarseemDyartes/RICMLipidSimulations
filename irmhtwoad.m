%% III: Head-Tail Simulation, No Adlayer
clear; clc;

%IRM Light Wavelength
lambda1 = 450e-9;
lambda2 = 532e-9; 
lambda3 = 610e-9;

%Refractive Indices
n0 = 1.52; %Substrate refractive index
n1 = 1.49; %Phosphocholine lipid head group refractive index
n2 = 1.44; %Lipid tail group refractive index (octadecane used in model)
n3 = 1.33; %PBS top solution refractive index

%Distances and Wavelength
m = 0:30;
npoint=max(size(m));
d1 = 0.5e-9; %Phosphocholine head thickness
d2 = 4e-9; %Tail thickness

%Fresnel Reflectance Coefficients
r01 = (n0-n1)/(n0+n1); %Substrate to head
r12 = (n1-n2)/(n1+n2); %Head to tail
r21 = -r12; %Tail to head
r13 = (n1-n3)/(n1+n3); %Head to top solution

%Fresnel Transmittance Coefficients
t01 = 2*n0/(n0+n1);  %Substrate to head
t12 = 2*n1/(n1+n2); %Head to tail
t21 = 2*n2/(n1+n2); %Tail to head
t13 = 2*n1/(n1+n3); %Head to top solution 

%Transfer Matrices between layers
M01 = [1,r01;r01,1]/t01; %Substrate to head
M12 = [1,r12;r12,1]/t12; %Head to tail
M21 = [1,r21;r21,1]/t21; %Tail to head
M13 = [1,r13;r13,1]/t13; %Head to top solution 

%Transfer Matrices within layers - 450 nm
phi1_l1 = 2*pi*n1*d1/lambda1; %Head
phi2_l1 = 2*pi*n2*d2/lambda1; %Tail
M1_l1 = [exp(1i*phi1_l1),0;0,exp(-1*1i*phi1_l1)]; %Head
M2_l1 = [exp(1i*phi2_l1),0;0,exp(-1*1i*phi2_l1)]; %Tail

%Transfer Matrices within layers - 532 nm
phi1_l2 = 2*pi*n1*d1/lambda2; %Head
phi2_l2 = 2*pi*n2*d2/lambda2; %Tail
M1_l2 = [exp(1i*phi1_l2),0;0,exp(-1*1i*phi1_l2)]; %Head
M2_l2 = [exp(1i*phi2_l2),0;0,exp(-1*1i*phi2_l2)]; %Tail

%Transfer Matrices within layers - 610 nm
phi1_l3 = 2*pi*n1*d1/lambda3; %Head
phi2_l3 = 2*pi*n2*d2/lambda3; %Tail
M1_l3 = [exp(1i*phi1_l3),0;0,exp(-1*1i*phi1_l3)]; %Head
M2_l3 = [exp(1i*phi2_l3),0;0,exp(-1*1i*phi2_l3)]; %Tail

%Coefficients
r=zeros(npoint,1);
t=r;
R=r;
T=r;
C=R;

%Light Pathway: m01 - j(m1 - m12 - m2 - m21 - m1) - m13
for j=1:npoint
	M_l1 = M01*[(M1_l1*M12*M2_l1*M21*M1_l1)^j]*M13;
	r_l1(j)=M_l1(2,1)/M_l1(1,1);
    t_l1(j)=1/M_l1(1,1);
	
	M_l2 = M01*[(M1_l2*M12*M2_l2*M21*M1_l2)^j]*M13;
	r_l2(j)=M_l2(2,1)/M_l2(1,1);
    t_l2(j)=1/M_l2(1,1);
	
	M_l3 = M01*[(M1_l3*M12*M2_l3*M21*M1_l3)^j]*M13;
	r_l3(j)=M_l3(2,1)/M_l3(1,1);
    t_l3(j)=1/M_l3(1,1);
end
for k=1:npoint
        R_l1(k)=norm(r_l1(k))^2;
        T_l1(k)=norm(t_l1(k))^2*n3/n0;
		
		R_l2(k)=norm(r_l2(k))^2;
        T_l2(k)=norm(t_l2(k))^2*n3/n0;
		
		R_l3(k)=norm(r_l3(k))^2;
        T_l3(k)=norm(t_l3(k))^2*n3/n0;
end

Rbg=(((n0-n3)./(n0+n3)).^2);
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
textct=sprintf('n_S_o_l=%4.3f\n n_H_e_a_d=%4.3f\n n_T_a_i_l%4.3f\n n_G_l_a_s_s=%4.3f\n',n3,n1,n2,n0);
% text('Position',[20 1 0]);
% text(min(xlim())+0.6*(max(xlim())-min(xlim())),0.8*min(ylim())+0.2*(max(ylim())-min(ylim())),textct);
text(7,0.9,textct);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
opfnhead=sprintf('nSol=%4.2f,nHead=%4.2f,nTail=%4.2f, wavelength=%4.1fnm, HeadTailNoAdlayer',n3,n1,n2,lambda1*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l1],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l1(2);
Mimg(412:612,412:612)=I_glass*CIRM_l1(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);

opfnhead=sprintf('nSol=%4.2f,nHead=%4.2f,nTail=%4.2f, wavelength=%4.1fnm, HeadTailNoAdlayer',n3,n1,n2,lambda2*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l2],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l2(2);
Mimg(412:612,412:612)=I_glass*CIRM_l2(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);

opfnhead=sprintf('nSol=%4.2f,nHead=%4.2f,nTail=%4.2f, wavelength=%4.1fnm, HeadTailNoAdlayer',n3,n1,n2,lambda3*1e9);
print([opfnhead,'.png'],'-dpng');
dlmwrite([opfnhead,'.dat'],[m,CIRM_l3],'\t');
I_glass=4000;
Mimg=I_glass*ones(1024,1024);
Mimg(212:812,212:812)=I_glass*CIRM_l3(2);
Mimg(412:612,412:612)=I_glass*CIRM_l3(3);
imwrite(uint16(Mimg),[opfnhead,'.tif']);