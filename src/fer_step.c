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

#include "fermi.h"

int ferstep_ST(void)
{

  if(fersrods(calcu.t)){
    PetscPrintf(FERMI_Comm,"Problem calculating control rods.\n");
    return 1;
  }
  if(ferass_ST()){
    PetscPrintf(FERMI_Comm,"Problem assembling the A & B.");
    return 1;
  }
  if(fersolv_ST()){
    PetscPrintf(FERMI_Comm,"Problem solving the steady state.\n");
    return 1;
  }

  /* Normalize the flux according to the power */
  if(fer_norm()){
    PetscPrintf(FERMI_Comm,"Problem normalizing the power.\n");
    return 1;
  }

  return 0;
}

int ferstep_TR(int step)
{

  if(fersrods(calcu.t)){
    PetscPrintf(FERMI_Comm,"Problem calculating control rods.\n");
    return 1;
  }
  if(ferass_TR(step)){
    PetscPrintf(FERMI_Comm,"Problem assembling the A & b.\n");
    return 1;
  }
  if(fersolv_TR()){
    PetscPrintf(FERMI_Comm,"Problem solving the transient steps.\n");
    return 1;
  }

  return 0;
}
