%this program perfoms and integral with gauss point evaluations

format long

xp=[ 0.211324865405187 0.211324865405187 0.211324865405187  ;
     0.211324865405187 0.788675134594813 0.211324865405187  ;
     0.788675134594813 0.788675134594813 0.211324865405187  ;
     0.788675134594813 0.211324865405187 0.211324865405187  ;
     0.211324865405187 0.211324865405187 0.788675134594813  ;
     0.211324865405187 0.788675134594813 0.788675134594813  ;
     0.788675134594813 0.788675134594813 0.788675134594813  ;
     0.788675134594813 0.211324865405187 0.788675134594813  ];
    
    
    
wp=[1 1 1 1 1 1 1 1]/8;

integral_1_1 = 0.0 ;
integral_1_2 = 0.0 ;
integral_1_3 = 0.0 ;
integral_1_4 = 0.0 ;
integral_1_5 = 0.0 ;
integral_1_6 = 0.0 ;
integral_1_7 = 0.0 ;
integral_1_8 = 0.0 ;

integral_2_1 = 0.0 ;
integral_2_2 = 0.0 ;
integral_2_3 = 0.0 ;
integral_2_4 = 0.0 ;
integral_2_5 = 0.0 ;
integral_2_6 = 0.0 ;
integral_2_7 = 0.0 ;
integral_2_8 = 0.0 ;

for i=1:8
    integral_1_1 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * wp(i);
    integral_1_2 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(1-xp(i,2))*(1-xp(i,3))) * wp(i);
    integral_1_3 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(xp(i,2))  *(1-xp(i,3))) * wp(i);
    integral_1_4 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(xp(i,2))  *(1-xp(i,3))) * wp(i);
    integral_1_5 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(1-xp(i,2))*(xp(i,3)))   * wp(i);
    integral_1_6 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(1-xp(i,2))*(xp(i,3)))   * wp(i);
    integral_1_7 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(xp(i,2))  *(xp(i,3)))   * wp(i);
    integral_1_8 += ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(xp(i,2))  *(xp(i,3)))   * wp(i);
    
    integral_2_1 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * wp(i);
    integral_2_2 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(1-xp(i,2))*(1-xp(i,3))) * wp(i);
    integral_2_3 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(xp(i,2))  *(1-xp(i,3))) * wp(i);
    integral_2_4 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(xp(i,2))  *(1-xp(i,3))) * wp(i);
    integral_2_5 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(1-xp(i,2))*(xp(i,3)))   * wp(i);
    integral_2_6 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(1-xp(i,2))*(xp(i,3)))   * wp(i);
    integral_2_7 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((xp(i,1))  *(xp(i,2))  *(xp(i,3)))   * wp(i);
    integral_2_8 += ((xp(i,1))*(1-xp(i,2))*(1-xp(i,3))) * ((1-xp(i,1))*(xp(i,2))  *(xp(i,3)))   * wp(i);
end

disp(integral_1_1)
disp(integral_1_2)
disp(integral_1_3)
disp(integral_1_4)
disp(integral_1_5)
disp(integral_1_6)
disp(integral_1_7)
disp(integral_1_8)
printf("\n")
disp(integral_2_1)
disp(integral_2_2)
disp(integral_2_3)
disp(integral_2_4)
disp(integral_2_5)
disp(integral_2_6)
disp(integral_2_7)
disp(integral_2_8)
printf("\n")

%evaluation of jacobian and dependance of the gauss point position 


x1=[ 1    1    1  ;
     3    1    1  ;  
     3    3    1  ;
     1    3    1  ;
     1    1    3  ;
     3    1    3  ;  
     3    3    3  ;
     1    3    3  ];
     
 
 % Each Jacobian will be evaluated at a diferent Gauss point
 jac1=zeros(3,3);
 jac2=zeros(3,3);
 jac3=zeros(3,3);
 jac4=zeros(3,3);
 jac5=zeros(3,3);
 jac6=zeros(3,3);
 jac7=zeros(3,3);
 jac8=zeros(3,3);
 jacx=zeros(3,3); 
 
 
 
%   grad_shape_3_8(2,1,0.5,0.5,0.5)
%   grad_shape_3_8(4,2,0.5,0.5,0.5)

 
 for i=1:3
    for j=1:3
       for n=1:8
          jac1(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(1,1),xp(1,2),xp(1,3)) ;
          jac2(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(2,1),xp(2,2),xp(2,3)) ;
          jac3(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(3,1),xp(3,2),xp(3,3)) ;
          jac4(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(4,1),xp(4,2),xp(4,3)) ;
          jac5(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(5,1),xp(5,2),xp(5,3)) ;
          jac6(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(6,1),xp(6,2),xp(6,3)) ;
          jac7(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(7,1),xp(7,2),xp(7,3)) ;
          jac8(i,j) += x1(n,j) * grad_shape_3_8(n,i,xp(8,1),xp(8,2),xp(8,3)) ;
          jacx(i,j) += x1(n,j) * grad_shape_3_8(n,i,0.5,0.5,0.5) ;
       end
    end
 end
 
 det(jac1)
 det(jac2)
 det(jac3)
 det(jac4)
 det(jac5)
 det(jac6)
 det(jac7)
 det(jac8)
 det(jacx)
 printf("\n")
 
% Triangle experiment
 
xp=zeros(4,3); 
wp=zeros(1,4);


 
xp=[0.1381966011250105   0.1381966011250105   0.1381966011250105 ;
    0.5854101966249685   0.1381966011250105   0.1381966011250105 ;
    0.1381966011250105   0.5854101966249685   0.1381966011250105 ;
    0.1381966011250105   0.1381966011250105   0.5854101966249685 ];
    
wp=[1 1 1 1]/(4*6);
 
x2=[ 0  5  5 ;
     5  0  5 ; 
     5  5  0 ;
     10 5  5 ];

printf("\n")
jac1=zeros(3,3);
jac2=zeros(3,3);
jac3=zeros(3,3);
jac4=zeros(3,3);
jacx=zeros(3,3);

 for i=1:3
   for j=1:3
      for n=1:4
          jac1(i,j) += x2(n,j) * grad_shape_3_4(n,i,xp(1,1),xp(1,2),xp(1,3)) ;
          jac2(i,j) += x2(n,j) * grad_shape_3_4(n,i,xp(2,1),xp(2,2),xp(2,3)) ;
          jac3(i,j) += x2(n,j) * grad_shape_3_4(n,i,xp(3,1),xp(3,2),xp(3,3)) ;
          jac4(i,j) += x2(n,j) * grad_shape_3_4(n,i,xp(4,1),xp(4,2),xp(4,3)) ;
          jacx(i,j) += x2(n,j) * grad_shape_3_4(n,i,0.5,0.5,0.5) ;
       end
    end
 end
 
 jac1
 
 det(jac1)
 det(jac2)
 det(jac3)
 det(jac4)
 det(jacx)
 
 integral_shapes=zeros(4,4); 
 for i=1:4
   for j=1:4
      for g=1:4
            integral_shapes(i,j) += shape_3_4(i,xp(g,1),xp(g,2),xp(g,3))* shape_3_4(j,xp(g,1),xp(g,2),xp(g,3))*wp(g);
      end   
   end
 end
 integral_shapes
 
 
 integral_dshapes=zeros(4,4); 
 for i=1:4
  for j=1:4
     for g=1:4
       for d=1:3 % the dot product
            integral_dshapes(i,j) += grad_shape_3_4(i,d,xp(g,1),xp(g,2),xp(g,3))*grad_shape_3_4(j,d,xp(g,1),xp(g,2),xp(g,3))*wp(g);
       end
     end   
   end
 end 
 integral_dshapes

