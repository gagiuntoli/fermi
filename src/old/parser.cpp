/*
 *  This source code is part of Fermi: a finite element code
 *  to solve the neutron diffusion problem for nuclear reactor
 *  designs.
 *
 *  Copyright (C) - 2019 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fermi.hpp"

#define NBUF 128

int parse_input(void) {
  if (access(inputfile, F_OK) == -1) {
    PetscPrintf(FERMI_Comm, "parser.c: file %s not found.\n", inputfile);
    return 1;
  }

  if (parse_mesh()) return 1;
  if (parse_mats()) return 2;
  if (parse_boun()) return 5;
  if (parse_outp()) return 7;
  return 0;
}

int parse_mesh(void) {
  FILE *file = fopen(inputfile, "r");
  char *data, buf[NBUF];
  int fl = 0, com = 0, ln = 0;

  while (fgets(buf, NBUF, file)) {
    ln++;
    data = strtok(buf, " \n");
    if (data) {
      if (!strcmp(data, "$Mesh")) {
        fl = 1;
      } else if (!strcmp(data, "$EndMesh")) {
        if (com == 7) {
          return 0;
        } else if (nproc == 1 && com >= 1) {
          return 0;
        } else if (nproc > 1 && (!(com & 4) || !(com & 8))) {
          PetscPrintf(FERMI_Comm, "parser.c:part file NF.\n");
          return 1;
        } else {
          PetscPrintf(FERMI_Comm, "parser.c:$Mesh sect BF.\n");
          return 1;
        }
      }
      if (fl == 2) {
        if (strcmp(data, "mesh_file") == 0) {
          data = strtok(NULL, " \n");
          if (!data) {
            PetscPrintf(FERMI_Comm, "parser.c:meshfile line %d.\n", ln);
            return 1;
          }
          strcpy(meshfile, data);
          if (access(meshfile, F_OK) == -1) {
            PetscPrintf(FERMI_Comm, "parser.c:%s NF.\n", meshfile);
            return 1;
          }
          com = com | 1;
        } else if (strcmp(data, "parfile_e") == 0) {
          data = strtok(NULL, " \n");
          if (!data) {
            PetscPrintf(FERMI_Comm, "parser.c:partfile line %d.\n", ln);
            return 1;
          }
          strcpy(epartfile, data);
          if (access(epartfile, F_OK) == -1) {
            PetscPrintf(FERMI_Comm, "parser.c:%s NF.\n", epartfile);
            return 1;
          }
          com = com | 2;
        } else if (strcmp(data, "parfile_n") == 0) {
          data = strtok(NULL, " \n");
          if (!data) {
            PetscPrintf(FERMI_Comm, "parser.c:partfile line %d.\n", ln);
            return 1;
          }
          strcpy(npartfile, data);
          if (access(npartfile, F_OK) == -1) {
            PetscPrintf(FERMI_Comm, "parser.c:%s NF.\n", npartfile);
            return 1;
          }
          com = com | 4;
        } else if (strcmp(data, "#") != 0) {
          PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
          return 1;
        }
      }
      if (fl == 1) fl = 2;
    }
  }
  fclose(file);
  return 1;
}

int parse_mats(void) {
  FILE *file = fopen(inputfile, "r");
  char *data, buf[NBUF], bufcpy[NBUF];
  int i, fl = 0, ln = 0, com = 0;
  pvl_t mat;

  while (fgets(buf, NBUF, file)) {
    ln++;
    strcpy(bufcpy, buf);
    data = strtok(buf, " \n");
    if (data) {
      if (!strcmp(data, "$Xs")) {
        fl = 1;
      } else if (!strcmp(data, "$EndXs")) {
        if (!list_mater.sizelist) {
          PetscPrintf(FERMI_Comm, "parser.c:no materials specified.\n");
          return 1;
        }
        return 0;
      }
      if (fl == 2) {
        if (!strcmp(data, "egn")) {
          data = strtok(NULL, " \n");
          if (!data) {
            PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
            return 1;
          }
          egn = atoi(data);
          if (egn < 1) {
            PetscPrintf(FERMI_Comm, "parser.c:egn should positive at line %d.\n", ln);
            return 1;
          }
          veloc = (double *)calloc(egn, sizeof(double));
          if (egn == 1) {
            nxs_mat = 5 * egn;
          } else {
            nxs_mat = 5 * egn + egn * (egn - 1);
          }
          com = com | 1;
        } else if (!strcmp(data, "pgn")) {
          data = strtok(NULL, " \n");
          if (!data) {
            PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
            return 1;
          }
          pgn = atoi(data);
          if (pgn < 1) {
            PetscPrintf(FERMI_Comm, "parser.c:pgn should positive at line %d.\n", ln);
            return 1;
          }
          beta = (double *)calloc(pgn, sizeof(double));
          lambda = (double *)calloc(pgn, sizeof(double));
          chi = (double *)calloc(pgn, sizeof(double));
          com = com | 2;
        } else if (!strcmp(data, "vel")) {
          for (i = 0; i < egn; i++) {
            data = strtok(NULL, " \n");
            if (!data) {
              PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
              return 1;
            }
            veloc[i] = atof(data);
          }
          com = com | 4;
        } else if (!strcmp(data, "kyn")) {
          for (i = 0; i < pgn; i++) {
            data = strtok(NULL, " \n");
            if (!data) {
              PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
              return 1;
            }
            beta[i] = atof(data);
          }
          for (i = 0; i < pgn; i++) {
            data = strtok(NULL, " \n");
            if (!data) {
              PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
              return 1;
            }
            lambda[i] = atof(data);
          }
          for (i = 0; i < pgn; i++) {
            data = strtok(NULL, " \n");
            if (!data) {
              PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
              return 1;
            }
            chi[i] = atof(data);
          }
          com = com | 8;
        } else if (data[0] != '#') {
          if ((com & 1) != 1) {
            PetscPrintf(FERMI_Comm, "parser.c:egn should be before xs values at line %d.\n", ln);
            return 1;
          }
          mat.D = (double *)calloc(egn, sizeof(double));
          mat.xs_a = (double *)calloc(egn, sizeof(double));
          mat.xs_s = (double *)calloc(egn * (egn - 1), sizeof(double));
          mat.nxs_f = (double *)calloc(egn, sizeof(double));
          mat.exs_f = (double *)calloc(egn, sizeof(double));
          mat.chi = (double *)calloc(egn, sizeof(double));
          if (parse_material(bufcpy, &mat)) {
            PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
            return 1;
          }
          if (list_insert_se(&list_mater, (void *)&mat)) return 1;
        }
      }
      if (fl == 1) fl = 2;
    }
  }
  fclose(file);
  return 1;
}

int parse_boun(void) {
  FILE *file = fopen(inputfile, "r");
  char *data, buf[NBUF], bufcpy[NBUF];
  int fl = 0, ln = 0;
  bound_t bou;

  while (fgets(buf, NBUF, file)) {
    ln++;
    strcpy(bufcpy, buf);
    data = strtok(buf, " \n");
    if (data) {
      if (!strcmp(data, "$Boundary")) {
        fl = 1;
      } else if (!strcmp(data, "$EndBoundary")) {
        if (!list_bound.sizelist) {
          PetscPrintf(FERMI_Comm, "parser.c: No boundaries in %s.\n", inputfile);
          return 1;
        }
        return 0;
      }
      if (fl == 2) {
        if (data[0] != '#') {
          if (parse_boundary(bufcpy, &bou)) {
            PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
            return 1;
          }
          list_insert_se(&list_bound, (void *)&bou);
        }
      }
      if (fl == 1) fl = 2;
    }
  }
  fclose(file);
  return 1;
}

int parse_material(char *buff, pvl_t *mat) {
  char *data;
  int i;

  if (!strlen(buff)) return 1;
  data = strtok(buff, " \n");
  if (!data) return 1;
  if (data[0] != '\"' || data[strlen(data) - 1] != '\"') return 1;
  strcpy(mat->name, data);

  /* second column is 0|1 which means if the mat has or not has precursors */
  data = strtok(NULL, " \n");
  if (!data) return 1;
  if (atoi(data) == 1 || atoi(data) == 0) {
    mat->hasprec = atoi(data);
  }

  /* diffusion coeficients */
  for (i = 0; i < egn; i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->D[i] = atof(data);
  }

  /* absortion cross sections */
  for (i = 0; i < egn; i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->xs_a[i] = atof(data);
  }

  /* scattering cross sections */
  for (i = 0; i < egn * (egn - 1); i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->xs_s[i] = atof(data);
  }

  /* fission cross sections */
  for (i = 0; i < egn; i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->nxs_f[i] = atof(data);
  }

  /* energy cross sections */
  for (i = 0; i < egn; i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->exs_f[i] = atof(data);
  }

  /* fission spectrum */
  for (i = 0; i < egn; i++) {
    data = strtok(NULL, " \n");
    if (!data) return 1;
    mat->chi[i] = atof(data);
  }
  return 0;
}

int parse_outp(void) {
  /*
     Parse for $Output word on parser
     it could be repeated, each one has a different kind
     which means what they want to output
   */

  FILE *file = fopen(inputfile, "r");
  char *data, buf[NBUF], buf_a[NBUF];
  int fl, ln, i;
  output_t output;

  ln = 0;
  fl = 0;
  while (fgets(buf, NBUF, file)) {
    ln++;
    strcpy(buf_a, buf);
    data = strtok(buf_a, " \n");

    if (data) {  // if the line is not empty

      if (!strcmp(data, "$Output")) {
        fl = 1;
      } else if (!strcmp(data, "$EndOutput")) {
        if (output.kind == 1) {
          list_insertlast(&list_outpu, (void *)&output);

        } else if (output.kind == 2) {
          list_insertlast(&list_outpu, (void *)&output);

        } else {
          PetscPrintf(FERMI_Comm, "parser.c:BF line %d.\n", ln);
          return 1;
        }
        fl = 0;
      }

      if (fl == 2 && data[0] != '#') {
        // we are in the line after $Output
        if (get_int(buf, "kind", &output.kind)) return 1;

        switch (output.kind) {
          case 1:
            break;

          case 2:

            /*
               kind = 2

               power on different physical entities on ASCII file
               file <file.dat>
               nphy <num_of_phys>
               "phys_1" "phys_2" ... "phys_n"

             */

            if (!fgets(buf, NBUF, file)) return 1;
            ln++;

            if (get_char(buf, "file", output.kind_2.file)) return 1;

            if (!fgets(buf, NBUF, file)) return 1;
            ln++;

            if (get_int(buf, "nphy", &output.kind_2.nphy)) return 1;
            output.kind_2.phys = (char **)malloc(output.kind_2.nphy * sizeof(char *));
            output.kind_2.ids = (int *)malloc(output.kind_2.nphy * sizeof(int *));
            output.kind_2.pow = (double *)malloc(output.kind_2.nphy * sizeof(double *));
            for (i = 0; i < output.kind_2.nphy; i++) {
              output.kind_2.phys[i] = (char *)malloc(16 * sizeof(char));
            }

            if (!fgets(buf, NBUF, file)) return 1;
            ln++;

            // now we read the physical entities names
            data = strtok(buf, " \n");
            i = 0;
            while (i < output.kind_2.nphy && data != NULL) {
              strcpy(output.kind_2.phys[i], data);
              data = strtok(NULL, " \n");
              i++;
            }
            if (i != output.kind_2.nphy) {
              return 1;
            }

            break;

          default:
            break;
        }
      }
      if (fl == 1) fl = 2;
    }
  }
  fclose(file);
  return 0;
}

int get_int(char *buf, const char *name, int *a) {
  /*
     Looks in "buf" for :
     "name" <(int) a>
     returns 1 if this is not correct
             0 if this is 0K
   */

  char *data;

  data = strtok(buf, " \n");

  if (!data) {
    PetscPrintf(FERMI_Comm, "%s expected.\n", name);
    return 1;
  }

  if (strcmp(data, name)) {
    PetscPrintf(FERMI_Comm, "%s expected.\n", name);
    return 1;
  }

  data = strtok(NULL, " \n");
  if (!data) {
    PetscPrintf(FERMI_Comm, "%s value expected.\n", name);
    return 1;
  }
  *a = atoi(data);

  return 0;
}

int get_char(char *buf, const char *name, char *a) {
  /*
     Looks in "buf" for :
     "name" <(char) a>
     returns 1 if this is not correct
             0 if this is 0K
   */

  char *data;

  data = strtok(buf, " \n");

  if (!data) {
    PetscPrintf(FERMI_Comm, "%s expected.\n", name);
    return 1;
  }

  if (strcmp(data, name)) {
    PetscPrintf(FERMI_Comm, "%s expected.\n", name);
    return 1;
  }

  data = strtok(NULL, " \n");
  if (!data) {
    PetscPrintf(FERMI_Comm, "%s value expected.\n", name);
    return 1;
  }
  strcpy(a, data);

  return 0;
}

int parse_boundary(char *buff, bound_t *bou) {
  char *data;

  if (!strlen(buff)) return 1;

  data = strtok(buff, " \n");
  if (!data) return 1;
  if (data[0] != '\"' || data[strlen(data) - 1] != '\"') return 1;
  strcpy(bou->name, data);

  data = strtok(NULL, " \n");
  if (!data) return 1;
  bou->order = atoi(data);

  data = strtok(NULL, " \n");
  if (!data) return 1;
  bou->kind = atoi(data);

  list_init(&bou->nodeL, sizeof(int), cmp_int);
  list_init(&bou->elemsL, sizeof(int), cmp_int);

  return 0;
}

int cmp_bou(void *a, void *b) {
  if (((bound_t *)a)->order > ((bound_t *)b)->order) {
    return 1;
  } else if (((bound_t *)a)->order == ((bound_t *)b)->order) {
    return 0;
  } else {
    return -1;
  }
}

int cmp_mat(void *a, void *b) { return (strcmp(((pvl_t *)a)->name, ((pvl_t *)b)->name) == 0) ? 0 : -1; }
