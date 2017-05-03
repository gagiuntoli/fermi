/* Material laws available
 * 
 * 
 */


#ifndef MATLAW_H
#define MATLAW_H

typedef int (*matlaw_t)(double **Cta, double *Stra, double *Stre, void *param, int dim, int fl, int *nparam);
int matlaw_ei(double **Cta, double *Stra, double *Stre, void *param, int dim, int fl, int *nparam);

#endif
