% code for testing elemental matrix building
% in 2D

D1=1.500; D2=0.400; XSA1=0.01012;  XSA2=0.080032;  NXSF1=0.000; NXSF2=0.135; CHI1=1.000;  CHI2=0.000;  XS12=0.000; XS21=0.020;  

n1=[130 0];
n2=[150 0];
n3=[139.95473251028901 9.8333333333104687];


jac=zeros(2,2);
sh=zeros(3,3);
dsh=zeros(3,2);
v1=zeros(2);
v2=zeros(2);

sh(1,1)=1-1/6-1/6; sh(1,2)=1-1/6-2/3; sh(1,3)=1-1/6-2/3; 
sh(2,1)=1/6;       sh(2,2)=2/3;       sh(2,3)=1/6;
sh(3,1)=1/6;       sh(3,2)=1/6;       sh(3,3)=2/3;

dsh(1,1)=-1.0; dsh(1,2)=-1.0; 
dsh(2,1)=+1.0; dsh(2,2)=+0.0; 
dsh(3,1)=+0.0; dsh(3,2)=+1.0; 

jac(1,1)=n1(1)*dsh(1,1)+n2(1)*dsh(2,1)+n3(1)*dsh(3,1); jac(1,2)=n1(2)*dsh(1,1)+n2(2)*dsh(2,1)+n3(2)*dsh(3,1);
jac(2,1)=n1(1)*dsh(1,2)+n2(1)*dsh(2,2)+n3(1)*dsh(3,2); jac(2,2)=n1(2)*dsh(1,2)+n2(2)*dsh(2,2)+n3(2)*dsh(3,2);

v1=inv(jac)*dsh(1,:)';
v2=inv(jac)*dsh(1,:)';

Ae=zeros(6,6);
Be=zeros(6,6);


for i=1:1
   Ae(1,1)+=(D1*v1'*v1 + (XSA1 + XS21)*sh(1,i)*sh(1,i) )*det(jac)/6; 
   Be(1,2)+= CHI1 * NXSF2*sh(1,i)*sh(1,i) *det(jac)/6;    Be(1,4)+= CHI1 * NXSF2*sh(2,i)*sh(1,i) *det(jac)/6; 

end