/*
 * elemental.c - RUTINES FOR ELEMENTAL MATRICES' BUILDING
 * 
 * Author : Guido Giuntoli
 *   
 */

#include "types.h"
#include "globals.h"

void elemental_matrix(int e, double *Ae, double *Be){
  
  int gp, i, j, g, m, g1, g2, d, mv, p;  
  double detj = 0.0;  
  double v1[MAX_DIM], v2[MAX_DIM], jac[MAX_DIM][MAX_DIM];
  
  shlist_t * shl;  
  
  shl = found_shape(mesh->elem[e].npe);
  
  for(i=0;i<mesh->elem[e].npe;i++){
    for(d=0;d<mesh->dim;d++){
      shl->coords[i][d]=mesh->elem[e].node[i].coord[d];
    }
  }
  
  if(egn == 1){
    
    for(i=0;i<mesh->elem[e].npe;i++){
      for(j=0;j<mesh->elem[e].npe;j++){
	Ae[ i*mesh->elem[e].npe*egn + j ]=0.0;
	Be[ i*mesh->elem[e].npe*egn + j ]=0.0;
	for(gp=0; gp<shl->gpn; gp++){
	  
	  detj = built_jacobian(shl->coords, shl->dsh, jac, mesh->dim, mesh->elem[e].npe, gp);
	  detj = fabs(detj);
	  memset(v1, 0.0, sizeof(double) * mesh->dim);
	  memset(v2, 0.0, sizeof(double) * mesh->dim);
	  
	  for(d=0;d<mesh->dim;d++){ 
	    for(mv=0;mv<mesh->dim;mv++){
	      v1[d]+= jac[d][mv]*shl->dsh[j % mesh->elem[e].npe][mv][gp];
	      v2[d]+= jac[d][mv]*shl->dsh[i % mesh->elem[e].npe][mv][gp];
	    }
	  }
	  
	  Ae[ i*mesh->elem[e].npe*egn + j ] += ( mesh->elem[e].D[0] * dot( v1, v2, mesh->dim) + mesh->elem[e].xs_a[0] * shl->sh[i % mesh->elem[e].npe][gp] * shl->sh[j % mesh->elem[e].npe][gp] ) * shl->wp[gp] * detj;
	  Be[ i*mesh->elem[e].npe*egn + j ] += mesh->elem[e].nxs_f[0] * shl->sh[j % mesh->elem[e].npe][gp] * shl->sh[i % mesh->elem[e].npe][gp] * shl->wp[gp] * detj;
	}
	
      }              
    }
    
  }else{
    
    m=0;
    p=0;
    for(i=0;i<mesh->elem[e].npe * egn;i++){
      if(i!=0 && (i % egn == 0)) m++;
      g=0;
      for(j=0;j<mesh->elem[e].npe * egn;j++){
	
	if(j!=0 && (j % egn == 0)) g++;
	Ae[ i*mesh->elem[e].npe*egn + j ]=0.0;
	Be[ i*mesh->elem[e].npe*egn + j ]=0.0;
	p=0;
	
	for(gp=0; gp<shl->gpn; gp++){
	  
	  detj = built_jacobian(shl->coords, shl->dsh, jac, mesh->dim, mesh->elem[e].npe, gp); 
	  detj=fabs(detj);
	  if( (j-i) % egn == 0 ){
	    
	    memset(v1, 0.0, sizeof(double) * mesh->dim);
	    memset(v2, 0.0, sizeof(double) * mesh->dim);
	    
	    for(d=0;d<mesh->dim;d++){ 
	      for(mv=0;mv<mesh->dim;mv++){
		v1[d]+= jac[d][mv]*shl->dsh[g][mv][gp];
		v2[d]+= jac[d][mv]*shl->dsh[m][mv][gp];
	      }
	    }
	    
	    Ae[ i*mesh->elem[e].npe*egn + j ] += ( mesh->elem[e].D[i % egn] * dot(v1,v2,mesh->dim)+ mesh->elem[e].xs_a[i % egn] * shl->sh[g][gp] * shl->sh[m][gp] ) * shl->wp[gp] * detj;
	    for(g1=0; g1<egn ; g1++){
	      for(g2=0; g2<egn ; g2++){
		if( (g1 != (i%egn)) && (g2 == (i%egn)) ){
		  if(g1>g2){
		    Ae[ i*mesh->elem[e].npe*egn + j ] += mesh->elem[e].xs_s[ g1 * (egn-1) + g2] * shl->sh[g][gp] * shl->sh[m][gp] * shl->wp[gp] * detj;
		  }else if(g1<g2){
		    Ae[ i*mesh->elem[e].npe*egn + j ] += mesh->elem[e].xs_s[ g1 * (egn-1) + g2 -1] * shl->sh[g][gp] * shl->sh[m][gp] * shl->wp[gp] * detj;
		  }
		}
	      }
	    }
	    
	  }else{   
	    
	    Ae[ i*mesh->elem[e].npe*egn + j ] -=  mesh->elem[e].xs_s[(i%egn)*(egn-1) + p ] * shl->sh[g][gp] * shl->sh[m][gp] * shl->wp[gp] * detj;
	    if(gp==(shl->gpn)-1)
	      p++; 
	    
	  }      
	  Be[ i*mesh->elem[e].npe*egn + j ] +=  mesh->elem[e].chi[i % egn] * mesh->elem[e].nxs_f[j % egn] * shl->sh[g][gp] * shl->sh[m][gp] * shl->wp[gp] * detj;
	}  
	
      } 
    }
    
  }  
  
}


void elemental_rhs(int e, double *be){

    // This rutine is for calculating the precursors contribution as a source term in the time-dependent equation
    
}


//void elemental_matrix(int e, double *Ae, double *Be){
//    
//    double detj;
//    int    n, i, j;
//    
//    /*we build the elemental matrix*/
//    if( mesh->dim == 1 ){
//        
//    }else if( mesh->dim == 2 ){
//        
//    }else if( mesh->dim == 3 ){ 
//        
//        double jac[3][3];
//        
//        if( mesh->elem[e].npe == 4 ){        
//            
//            /*******************************************************/
//            // TETRAHEDRON
//            /*******************************************************/
//            int    d, d1, d2;
//            double intsh[4][4], intdsh[4][4], v[4][3];
//            
//            /*******************************************************/
//            /* BUILD JACOBIAN MATRIX                               */
//            for(i=0;i<3;i++){
//                for(j=0;j<3;j++){
//                    jac[i][j]=0.0;
//                    for(n=0;n<mesh->elem[e].npe;n++){
//                        jac[i][j] += mesh->elem[e].node[n].coord[j] * dshape_3_4[n][i];
//                    }
//                }    
//            }
//            /*******************************************************/
//            detj = invmat_3( jac );
//            
//            for(d=0;d<4;d++){
//                for(d1=0;d1<3;d1++){
//                    v[d][d1] = 0.0;
//                    for(d2=0;d2<3;d2++){
//                        v[d][d1] +=  jac[d1][d2] * dshape_3_4[d][d2];
//                    }
//                } 
//            }
//            for(i=0;i<4;i++){
//                for(j=0;j<4;j++){
//                    if( i==j ){
//                        intsh[i][j] = detj/60;
//                    }else{
//                        intsh[i][j] = detj/120;
//                    }
//                    intdsh[i][j]=0.0;
//                    for(d=0;d<3;d++){
//                        intdsh[i][j] += v[i][d] * v[j][d];
//                    }
//                    intdsh[i][j] *= detj;
//                }                    
//            }                
//            /*****************************************************/  
//            
//            if( egn == 1 ){
//                
//                Ae[0]  = mesh->elem[e].D[0]*intdsh[0][0] + mesh->elem[e].xs_a[0]*intsh[0][0] ; Ae[1]  = mesh->elem[e].D[0]*intdsh[1][0] + mesh->elem[e].xs_a[0]*intsh[1][0] ; Ae[2]  = mesh->elem[e].D[0]*intdsh[2][0] + mesh->elem[e].xs_a[0]*intsh[2][0] ; Ae[3]  = mesh->elem[e].D[0]*intdsh[3][0] + mesh->elem[e].xs_a[0]*intsh[3][0] ; 
//                Ae[4]  = mesh->elem[e].D[0]*intdsh[0][1] + mesh->elem[e].xs_a[0]*intsh[0][1] ; Ae[5]  = mesh->elem[e].D[0]*intdsh[1][1] + mesh->elem[e].xs_a[0]*intsh[1][1] ; Ae[6]  = mesh->elem[e].D[0]*intdsh[2][1] + mesh->elem[e].xs_a[0]*intsh[2][1] ; Ae[7]  = mesh->elem[e].D[0]*intdsh[3][1] + mesh->elem[e].xs_a[0]*intsh[3][1] ;
//                Ae[8]  = mesh->elem[e].D[0]*intdsh[0][2] + mesh->elem[e].xs_a[0]*intsh[0][2] ; Ae[9]  = mesh->elem[e].D[0]*intdsh[1][2] + mesh->elem[e].xs_a[0]*intsh[1][2] ; Ae[10] = mesh->elem[e].D[0]*intdsh[2][2] + mesh->elem[e].xs_a[0]*intsh[2][2] ; Ae[11] = mesh->elem[e].D[0]*intdsh[3][2] + mesh->elem[e].xs_a[0]*intsh[3][2] ;
//                Ae[12] = mesh->elem[e].D[0]*intdsh[0][3] + mesh->elem[e].xs_a[0]*intsh[0][3] ; Ae[13] = mesh->elem[e].D[0]*intdsh[1][3] + mesh->elem[e].xs_a[0]*intsh[1][3] ; Ae[14] = mesh->elem[e].D[0]*intdsh[2][3] + mesh->elem[e].xs_a[0]*intsh[2][3] ; Ae[15] = mesh->elem[e].D[0]*intdsh[3][3] + mesh->elem[e].xs_a[0]*intsh[3][3] ;
//                
//                Be[0]  = mesh->elem[e].nxs_f[0]*intsh[0][0] ;Be[1]  = mesh->elem[e].nxs_f[0]*intsh[1][0] ; Be[2]  = mesh->elem[e].nxs_f[0]*intsh[2][0] ;Be[3] = mesh->elem[e].nxs_f[0]*intsh[3][0] ;   
//                Be[4]  = mesh->elem[e].nxs_f[0]*intsh[0][1] ;Be[5]  = mesh->elem[e].nxs_f[0]*intsh[1][1] ; Be[6]  = mesh->elem[e].nxs_f[0]*intsh[2][1] ;Be[3] = mesh->elem[e].nxs_f[0]*intsh[3][1] ;
//                Be[8]  = mesh->elem[e].nxs_f[0]*intsh[0][2] ;Be[9]  = mesh->elem[e].nxs_f[0]*intsh[1][2] ; Be[10] = mesh->elem[e].nxs_f[0]*intsh[2][2] ;Be[3] = mesh->elem[e].nxs_f[0]*intsh[3][2] ;
//                Be[12] = mesh->elem[e].nxs_f[0]*intsh[0][3] ;Be[13] = mesh->elem[e].nxs_f[0]*intsh[1][3] ; Be[14] = mesh->elem[e].nxs_f[0]*intsh[2][3] ;Be[3] = mesh->elem[e].nxs_f[0]*intsh[3][3] ;                
//                
//            }else if( egn == 2 ){
//                
//                Ae[0]  =  mesh->elem[e].D[0] * intdsh[0][0] + mesh->elem[e].xs_r[0] * intsh[0][0] ;  Ae[1]  = -mesh->elem[e].xs_s[1] * intsh[0][0]                                       ;  /* */  Ae[2]  =  mesh->elem[e].D[0]   *intdsh[1][0] + mesh->elem[e].xs_r[0] * intsh[1][0] ; Ae[3]  = -mesh->elem[e].xs_s[1]*intsh[1][0] ;                                       /* */ Ae[4]  =  mesh->elem[e].D[0]  *intdsh[2][0] + mesh->elem[e].xs_r[0] * intsh[2][0] ;  Ae[5]  = -mesh->elem[e].xs_s[1]*intsh[2][0] ;                                       /* */ Ae[6]  =  mesh->elem[e].D[0]  *intdsh[3][0] + mesh->elem[e].xs_r[0] * intsh[3][0] ;  Ae[7]  = -mesh->elem[e].xs_s[1]*intsh[3][0] ;
//                Ae[8]  = -mesh->elem[e].xs_s[0]*intsh[0][0]                                       ;  Ae[9]  =  mesh->elem[e].D[1]   * intdsh[0][0] + mesh->elem[e].xs_r[1] * intsh[0][0] ;  /* */  Ae[10] = -mesh->elem[e].xs_s[0]* intsh[1][0]                                       ; Ae[11] =  mesh->elem[e].D[1] * intdsh[1][0] + mesh->elem[e].xs_r[1] * intsh[1][0];  /* */ Ae[12] = -mesh->elem[e].xs_s[0]*intsh[2][0] ;                                        Ae[13] =  mesh->elem[e].D[1]  *intdsh[2][0] + mesh->elem[e].xs_r[1] * intsh[2][0];  /* */ Ae[14] = -mesh->elem[e].xs_s[0]*intsh[3][0]                                       ;  Ae[15] =  mesh->elem[e].D[1]  *intdsh[3][0] + mesh->elem[e].xs_r[1] * intsh[3][0];    
//                Ae[16] =  mesh->elem[e].D[0] * intdsh[0][1] + mesh->elem[e].xs_r[0] * intsh[0][1] ;  Ae[17] = -mesh->elem[e].xs_s[1] * intsh[0][1]                                       ;  /* */  Ae[18] =  mesh->elem[e].D[0]   *intdsh[1][1] + mesh->elem[e].xs_r[0] * intsh[1][1] ; Ae[19] = -mesh->elem[e].xs_s[1]*intsh[1][1] ;                                       /* */ Ae[20] =  mesh->elem[e].D[0]  *intdsh[2][1] + mesh->elem[e].xs_r[0] * intsh[2][1] ;  Ae[21] = -mesh->elem[e].xs_s[1]*intsh[2][1] ;                                       /* */ Ae[22] =  mesh->elem[e].D[0]  *intdsh[3][1] + mesh->elem[e].xs_r[0] * intsh[3][1] ;  Ae[23] = -mesh->elem[e].xs_s[1]*intsh[3][1] ;
//                Ae[24] = -mesh->elem[e].xs_s[0]*intsh[0][1]                                       ;  Ae[25] =  mesh->elem[e].D[1]   * intdsh[0][1] + mesh->elem[e].xs_r[1] * intsh[0][1] ;  /* */  Ae[26] = -mesh->elem[e].xs_s[0]* intsh[1][1]                                       ; Ae[27] =  mesh->elem[e].D[1] * intdsh[1][1] + mesh->elem[e].xs_r[1] * intsh[1][1];  /* */ Ae[28] = -mesh->elem[e].xs_s[0]*intsh[2][1] ;                                        Ae[29] =  mesh->elem[e].D[1]  *intdsh[2][1] + mesh->elem[e].xs_r[1] * intsh[2][1];  /* */ Ae[30] = -mesh->elem[e].xs_s[0]*intsh[3][1]                                       ;  Ae[31] =  mesh->elem[e].D[1]  *intdsh[3][1] + mesh->elem[e].xs_r[1] * intsh[3][1];              
//                Ae[32] =  mesh->elem[e].D[0] * intdsh[0][2] + mesh->elem[e].xs_r[0] * intsh[0][2] ;  Ae[33] = -mesh->elem[e].xs_s[1] * intsh[0][2]                                       ;  /* */  Ae[34] =  mesh->elem[e].D[0]   *intdsh[1][2] + mesh->elem[e].xs_r[0] * intsh[1][2] ; Ae[35] = -mesh->elem[e].xs_s[1]*intsh[1][2] ;                                       /* */ Ae[36] =  mesh->elem[e].D[0]  *intdsh[2][2] + mesh->elem[e].xs_r[0] * intsh[2][2] ;  Ae[37] = -mesh->elem[e].xs_s[1]*intsh[2][2] ;                                       /* */ Ae[38] =  mesh->elem[e].D[0]  *intdsh[3][2] + mesh->elem[e].xs_r[0] * intsh[3][2] ;  Ae[39] = -mesh->elem[e].xs_s[1]*intsh[3][2] ;
//                Ae[40] = -mesh->elem[e].xs_s[0]*intsh[0][2]                                       ;  Ae[41] =  mesh->elem[e].D[1]   * intdsh[0][2] + mesh->elem[e].xs_r[1] * intsh[0][2] ;  /* */  Ae[42] = -mesh->elem[e].xs_s[0]* intsh[1][2]                                       ; Ae[43] =  mesh->elem[e].D[1] * intdsh[1][2] + mesh->elem[e].xs_r[1] * intsh[1][2];  /* */ Ae[44] = -mesh->elem[e].xs_s[0]*intsh[2][2] ;                                        Ae[45] =  mesh->elem[e].D[1]  *intdsh[2][2] + mesh->elem[e].xs_r[1] * intsh[2][2];  /* */ Ae[46] = -mesh->elem[e].xs_s[0]*intsh[3][2]                                       ;  Ae[47] =  mesh->elem[e].D[1]  *intdsh[3][2] + mesh->elem[e].xs_r[1] * intsh[3][2];   
//                Ae[48] =  mesh->elem[e].D[0] * intdsh[0][3] + mesh->elem[e].xs_r[0] * intsh[0][3] ;  Ae[49] = -mesh->elem[e].xs_s[1] * intsh[0][3]                                       ;  /* */  Ae[50] =  mesh->elem[e].D[0]   *intdsh[1][3] + mesh->elem[e].xs_r[0] * intsh[1][3] ; Ae[51] = -mesh->elem[e].xs_s[1]*intsh[1][3] ;                                       /* */ Ae[52] =  mesh->elem[e].D[0] * intdsh[2][3] + mesh->elem[e].xs_r[0] * intsh[2][3] ;  Ae[53] = -mesh->elem[e].xs_s[1]*intsh[2][3] ;                                       /* */ Ae[54] =  mesh->elem[e].D[0]  *intdsh[3][3] + mesh->elem[e].xs_r[0] * intsh[3][3] ;  Ae[55] = -mesh->elem[e].xs_s[1]*intsh[3][3] ;
//                Ae[56] = -mesh->elem[e].xs_s[0]*intsh[0][3]                                       ;  Ae[57] =  mesh->elem[e].D[1]   * intdsh[0][3] + mesh->elem[e].xs_r[1] * intsh[0][3] ;  /* */  Ae[58] = -mesh->elem[e].xs_s[0]* intsh[1][3]                                       ; Ae[59] =  mesh->elem[e].D[1] * intdsh[1][3] + mesh->elem[e].xs_r[1] * intsh[1][3];  /* */ Ae[60] = -mesh->elem[e].xs_s[0]*intsh[2][3] ;                                        Ae[61] =  mesh->elem[e].D[1]  *intdsh[2][3] + mesh->elem[e].xs_r[1] * intsh[2][3];  /* */ Ae[62] = -mesh->elem[e].xs_s[0]*intsh[3][3]                                       ;  Ae[63] =  mesh->elem[e].D[1]  *intdsh[3][3] + mesh->elem[e].xs_r[1] * intsh[3][3]; 
//                
//                Be[0]   =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[0][0]            ;  Be[1]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[0][0]                ;  /* */  Be[2]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[1][0]               ;  Be[3]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[1][0]            ;  /* */ Be[4]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[2][0]              ;  Be[5]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[2][0]             ;  /* */ Be[6]  = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[3][0]              ;  Be[7]  =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[3][0] ;
//                Be[8]   =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[0][0]            ;  Be[9]  = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[0][0]                ;  /* */  Be[10] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[1][0]               ;  Be[11] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[1][0]            ;  /* */ Be[12] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[2][0]              ;  Be[13] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[2][0]             ;  /* */ Be[14] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[3][0]              ;  Be[15] =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[3][0] ;
//                Be[16]  =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[0][1]            ;  Be[17] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[0][1]                ;  /* */  Be[18] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[1][1]               ;  Be[19] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[1][1]            ;  /* */ Be[20] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[2][1]              ;  Be[21] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[2][1]             ;  /* */ Be[22] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[3][1]              ;  Be[23] =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[3][1] ;
//                Be[24]  =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[0][1]            ;  Be[25] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[0][1]                ;  /* */  Be[26] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[1][1]               ;  Be[27] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[1][1]            ;  /* */ Be[28] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[2][1]              ;  Be[29] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[2][1]             ;  /* */ Be[30] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[3][1]              ;  Be[31] =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[3][1] ;
//                Be[32]  =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[0][2]            ;  Be[33] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[0][2]                ;  /* */  Be[34] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[1][2]               ;  Be[35] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[1][2]            ;  /* */ Be[36] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[2][2]              ;  Be[37] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[2][2]             ;  /* */ Be[38] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[3][2]              ;  Be[39] =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[3][2] ;
//                Be[40]  =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[0][2]            ;  Be[41] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[0][2]                ;  /* */  Be[42] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[1][2]               ;  Be[43] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[1][2]            ;  /* */ Be[44] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[2][2]              ;  Be[45] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[2][2]             ;  /* */ Be[46] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[3][2]              ;  Be[47] =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[3][2] ;
//                Be[48]  =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[0][3]            ;  Be[49] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[0][3]                ;  /* */  Be[50] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[1][3]               ;  Be[51] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[1][3]            ;  /* */ Be[52] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[2][3]              ;  Be[53] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[2][3]             ;  /* */ Be[54] = mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[0] * intsh[3][3]              ;  Be[55] =  mesh->elem[e].chi[0] * mesh->elem[e].nxs_f[1] * intsh[3][3] ;
//                Be[56]  =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[0][3]            ;  Be[57] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[0][3]                ;  /* */  Be[58] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[1][3]               ;  Be[59] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[1][3]            ;  /* */ Be[60] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[2][3]              ;  Be[61] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[2][3]             ;  /* */ Be[62] = mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[0] * intsh[3][3]              ;  Be[63] =  mesh->elem[e].chi[1] * mesh->elem[e].nxs_f[1] * intsh[3][3] ;
//                
//            } 
//            
//        } else if( mesh->elem[e].npe == 8 ){ 
//            /*******************************************************/
//            // HEXAHEDRON
//            /*******************************************************/
//            
//        } 
//        
//    }else{
//        ierr = PetscPrintf( PETSC_COMM_WORLD, "error in elemdiff.c with dimension dim = %d \n", mesh->dim );
//        exit(1);
//    }   
//    return;
//}