/*****************************************************************************
 *
 *  psi_init.c
 *
 *  Various initial states for electrokinetics.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2023 The University of Edinburgh
 *
 *  Contributing authors:
 *  Oliver Henrich (ohenrich@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
 
#include "pe.h"
#include "coords.h"
#include "psi_init.h"

int *fixed_potential = NULL;  
int line_count = 0;

char **potential_evolution = NULL;
int line_count_evolution = 0;

double POTENTIAL = 25.0; //0.00001665;

/*****************************************************************************
 *
 *  psi_init_uniform
 *
 *  Set the charge density for all species to be rho_el everywhere.
 *  The potential is initialised to zero.
 *
 *****************************************************************************/

int psi_init_uniform(psi_t * obj, double rho_el) {

  int ic, jc, kc, index;
  int nlocal[3];
  int n, nk;
  
  assert(obj);
  assert(rho_el >= 0.0);

  cs_nlocal(obj->cs, nlocal);
  psi_nk(obj, &nk);

  int ntotal[3];
  cs_ntotal(obj->cs, ntotal);
  fixed_potential = (int *) calloc(ntotal[X] * ntotal[Y] * ntotal[Z], sizeof(int));

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);

	psi_psi_set(obj, index, 0.0);

	for (n = 0; n < nk; n++) {
	  psi_rho_set(obj, index, n, rho_el);
	}

      }
    }
  }
  
  return 0;
}

/*****************************************************************************
 *
 * psi_init_gouy_chapman
 *
 *  Set rho(x = 1)  = + (1/2NyNz)
 *      rho(x = Lx) = + (1/2NyNz)
 *      rho         = - 1/(NyNz*(Nx-2)) + electrolyte
 *
 *  This sets up the system for Gouy-Chapman.
 *
 *  rho_el is the electrolyte (background) charge density.
 *  sigma is the sufrace charge density at the wall.
 *
 *****************************************************************************/

int psi_init_gouy_chapman(psi_t * obj, map_t * map, double rho_el,
			      double sigma) {

  int ic, jc, kc, index;
  int nlocal[3];
  int mpi_cartsz[3];
  int mpi_cartcoords[3];
  double rho_w, rho_i;
  double ltot[3];

  assert(obj);
  assert(map);

  cs_nlocal(obj->cs, nlocal);
  cs_ltot(obj->cs, ltot);
  cs_cartsz(obj->cs, mpi_cartsz);
  cs_cart_coords(obj->cs, mpi_cartcoords);

  int ntotal[3];
  cs_ntotal(obj->cs, ntotal);
  fixed_potential = (int *) calloc(ntotal[X] * ntotal[Y] * ntotal[Z], sizeof(int));

  /* wall surface charge density */
  rho_w = sigma;

  /* counter charge density */
  rho_i = rho_w * 2.0 *ltot[Y]*ltot[Z] / (ltot[Y]*ltot[Z]*(ltot[X] - 2.0));

  /* apply counter charges & electrolyte */
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);

	psi_psi_set(obj, index, 0.0);
	psi_rho_set(obj, index, 0, rho_el);
	psi_rho_set(obj, index, 1, rho_el + rho_i);

      }
    }
  }

  /* apply wall charges */
  if (mpi_cartcoords[X] == 0) {
    ic = 1;
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);
	map_status_set(map, index, MAP_BOUNDARY);

	psi_rho_set(obj, index, 0, rho_w);
	psi_rho_set(obj, index, 1, 0.0);

      }
    }
  }

  if (mpi_cartcoords[X] == mpi_cartsz[X] - 1) {
    ic = nlocal[X];
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);
	map_status_set(map, index, MAP_BOUNDARY);

	psi_rho_set(obj, index, 0, rho_w);
	psi_rho_set(obj, index, 1, 0.0);

      }
    }
  }

  map_halo(map);
  map_pm_set(map, 1);

  return 0;
}

/*****************************************************************************
 *
 *  psi_init_liquid_junction
 *
 *  This follows Mafe et al., so we set for two species:
 *
 *    rho_left  = rho_el + delta_el / 2
 *    rho_right = rho_el - delta_el / 2
 *
 *  where left and right are separated by the half way point in the
 *  x-direction.
 *
 *  We should have delta_el << rho_el.
 *
 *****************************************************************************/

int psi_init_liquid_junction(psi_t * obj, double rho_el, double delta_el) {

  int ic, jc, kc, index;
  int ntotal[3];
  int nlocal[3], noff[3];

  assert(obj);

  cs_nlocal(obj->cs, nlocal);
  cs_ntotal(obj->cs, ntotal);
  cs_nlocal_offset(obj->cs, noff);

  fixed_potential = (int *) calloc(ntotal[X] * ntotal[Y] * ntotal[Z], sizeof(int));


  /* Set electrolyte densities */

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);

	psi_psi_set(obj, index, 0.0);

	if ((1 <= noff[X] + ic) && (noff[X] + ic <= ntotal[X]/2)) {
	  psi_rho_set(obj, index, 0, rho_el + 0.5*delta_el);
	  psi_rho_set(obj, index, 1, rho_el + 0.5*delta_el);
	}
	else {
	  psi_rho_set(obj, index, 0, rho_el - 0.5*delta_el);
	  psi_rho_set(obj, index, 1, rho_el - 0.5*delta_el);
	}
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 * psi_init_sigma
 *
 *  This sets up surface charges as specified in the porous media file.
 *
 *****************************************************************************/

int psi_init_sigma(psi_t * psi, map_t * map) {

  int ic, jc, kc, index;
  int nlocal[3];
  int ntotal[3];
  double sigma; /* point charge or surface charge density */

  
  assert(psi);
  assert(map);

  cs_nlocal(psi->cs, nlocal);
  cs_ntotal(psi->cs, ntotal);

  fixed_potential = (int *) calloc(ntotal[X] * ntotal[Y] * ntotal[Z], sizeof(int));

  //for (int i = 0; i < ntotal[X] * ntotal[Y] * ntotal[Z]; i++) {
    //printf("initialised fixed_potential[%d] = %d \n", i, fixed_potential[i]);
  //}

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(psi->cs, ic, jc, kc);
	map_data(map, index, &sigma);

	psi_psi_set(psi, index, 0.0);

	if (sigma) {
	  if (sigma > 0) {
	    psi_rho_set(psi, index, 0, sigma);
	    psi_rho_set(psi, index, 1, 0);
	  }
	  if (sigma < 0) {
	    psi_rho_set(psi, index, 0, 0);
	    psi_rho_set(psi, index, 1, -sigma);
	  }
	}

      }
    }
  }

  map_halo(map);

  return 0;
}

int psi_init_sigma_fixed_potential(psi_t * psi, map_t * map) {

  int ic, jc, kc, index;
  int mpi_cartsz[3];
  int mpicoords[3];
  int ntotal[3];

  printf("Starting psi_init_sigma_fixed_potential \n");

  cs_cartsz(psi->cs, mpi_cartsz);
  cs_cart_coords(psi->cs, mpicoords);
  cs_ntotal(psi->cs, ntotal);
  //int *fixed_potential;
  //int line_count;
  int counter = mpicoords[X] * (ntotal[Y] * ntotal[Z]) * ntotal[X] / mpi_cartsz[X] + mpicoords[Y] * ntotal[Z] * ntotal[Y] / mpi_cartsz[Y] + mpicoords[Z] * ntotal[Z] / mpi_cartsz[Z];
  int nlocal[3];
  int shift_z = ntotal[Z] * (mpi_cartsz[Z] - 1) / mpi_cartsz[Z];
  int shift_y = ntotal[Y] * ntotal[Z] * (mpi_cartsz[Y] - 1) / mpi_cartsz[Y];

  double sigma; /* point charge or surface charge density */
  fixed_potential = read_map_file("potential_map.txt", &line_count);
  //double * __restrict__ psidata = psi->psi->data;
  
  assert(psi);
  assert(map);

  cs_nlocal(psi->cs, nlocal);
  
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	      index = cs_index(psi->cs, ic, jc, kc);
	      map_data(map, index, &sigma);

        if (counter >= line_count) {
          printf("Errore: counter %d supera line_count %d\n", counter, line_count);
          exit(1);
          }
        
        if (fixed_potential[counter] == 1){
          psi->psi->data[addr_rank0(psi->nsites, index)] = POTENTIAL;
          counter++;  
        } 
        else if (fixed_potential[counter] == -1){
          psi->psi->data[addr_rank0(psi->nsites, index)] = -POTENTIAL;
          counter++; 
        }
	      else if (fixed_potential[counter] == 2){
          psi->psi->data[addr_rank0(psi->nsites, index)] = 0.0;
          counter++; 
        }
        else if (fixed_potential[counter] == 0){
	      psi->psi->data[addr_rank0(psi->nsites, index)] = 0.0; //commenta se vuoi farlo partire da 0
          counter++;
        }
        else {
          printf("ERROR READING potential_map \n");
	        exit(0);
          //psi->psi->data[addr_rank0(psi->nsites, index)] = 0.0; //commenta se vuoi farlo partire da 0
          //counter++;
        }
        //printf("row: %d psi: %f \n", ic + jc + kc, psi->psi->data[addr_rank0(psi->nsites, index)]);
	      //psi_psi_set(psi, index, 1.0);

	      if (sigma) {
	        if (sigma > 0) {
	          psi_rho_set(psi, index, 0, sigma);
	          psi_rho_set(psi, index, 1, 0);
	       }
	       if (sigma < 0) {
	          psi_rho_set(psi, index, 0, 0);
	          psi_rho_set(psi, index, 1, -sigma);
	        }
	      }

      }
      counter = counter + shift_z;
    }
    counter = counter + shift_y;
  }

  map_halo(map);
  

  return 0;
}



//////////////////////////////////////////
//
//
//  funzione per copiare i valori fissi del potenziale
//
//
//////////////////////////////////////////

int* read_map_file(const char *filename, int *line_count) {
  FILE *file;
  char line[1024];
  int *array;
  int value;

  // Apre il file per la lettura
  file = fopen(filename, "r");
  if (file == NULL) {
      perror("Errore nell'aprire il file");
      return NULL;
  }

  // Conta il numero di righe nel file
  *line_count = 0;
  while (fgets(line, sizeof(line), file) != NULL) {
      (*line_count)++;
  }

  // Riavvolge il puntatore del file all'inizio
  rewind(file);

  // Alloca l'array per contenere i valori (0 o 1)
  array = (int *)malloc(*line_count * sizeof(int));
  if (array == NULL) {
      perror("Errore nell'allocare la memoria");
      fclose(file);
      return NULL;
  }

  // Legge il file linea per linea e salva i valori nell'array
  int i = 0;
  while (fgets(line, sizeof(line), file) != NULL) {
      // Converte il valore della riga in un intero (0 o 1)
      if (sscanf(line, "%d", &value) == 1){ //&& (value == 0 || value == 1 || value == -1 || value == 2)) {
          array[i] = value;
          i++;
      } else {
          printf("Errore nel leggere il valore dalla riga: %s", line);
      }
  }

  // Chiude il file
  fclose(file);

  return array;
}

//////////////////////////////////////////
//
//
//  funzione per registrare i valori in cui evolve il potenziale
//
//
//////////////////////////////////////////

char **read_evolution_file(const char *filename, int *line_count) {
    FILE *file;
    char line[1024];
    char **array;

    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Errore nell'aprire il file");
        return NULL;
    }

    *line_count = 0;
    while (fgets(line, sizeof(line), file) != NULL) {
        (*line_count)++;
    }

    rewind(file);

    array = (char **)malloc(*line_count * sizeof(char *));
    if (array == NULL) {
        perror("Errore nell'allocare la memoria");
        fclose(file);
        return NULL;
    }

    int i = 0;
    while (fgets(line, sizeof(line), file) != NULL && i < *line_count) {
        line[strcspn(line, "\r\n")] = '\0';  // rimuove newline e carriage return
        if (strcmp(line, "t_0+") == 0 || strcmp(line, "t_1+") == 0 ||
            strcmp(line, "t_2+") == 0 || strcmp(line, "t_3+") == 0 ||
            strcmp(line, "t_0-") == 0 || strcmp(line, "t_1-") == 0 ||
            strcmp(line, "t_2-") == 0 || strcmp(line, "t_3-") == 0 || strcmp(line, "0") == 0) {
            array[i] = strdup(line);  // copia la stringa
            i++;
        } else {
            printf("Errore nel leggere il valore dalla riga: %s\n", line);
        }
    }

    /*
      array[i] = strdup(line); // copia la stringa
      array[i][strcspn(array[i], "\n")] = 0; // rimuove newline
      i++;
    */

    fclose(file);
    *line_count = i;  // aggiornamento reale in caso di linee scartate
    return array;
}


//////////////////////////////////////////
//
//
//  funzione per evolvere il potenziale
//
//
//////////////////////////////////////////

int psi_evolve_potential(psi_t * psi, map_t * map, int contatore_update_potenziale) {

  int ic, jc, kc, index;

  int nlocal[3];
  int ntotal[3];
  int mpi_cartsz[3];
  int mpicoords[3];

  cs_nlocal(psi->cs, nlocal);
  cs_ntotal(psi->cs, ntotal);
  cs_cartsz(psi->cs, mpi_cartsz);
  cs_cart_coords(psi->cs, mpicoords);

  int shift_z = ntotal[Z] * (mpi_cartsz[Z] - 1) / mpi_cartsz[Z];
  int shift_y = ntotal[Y] * ntotal[Z] * (mpi_cartsz[Y] - 1) / mpi_cartsz[Y];

  //int *fixed_potential;
  //int line_count;
  
  int counter = mpicoords[X] * (ntotal[Y] * ntotal[Z]) * ntotal[X] / mpi_cartsz[X] + mpicoords[Y] * ntotal[Z] * ntotal[Y] / mpi_cartsz[Y] + mpicoords[Z] * ntotal[Z] / mpi_cartsz[Z];
  double sigma; // point charge or surface charge density 
  potential_evolution = read_evolution_file("potential_evolution.txt", &line_count);
  //double * __restrict__ psidata = psi->psi->data;
  double potential_zero = 0.0;
  assert(psi);
  assert(map);

  cs_nlocal(psi->cs, nlocal);
  
  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	      index = cs_index(psi->cs, ic, jc, kc);
        
	      map_data(map, index, &sigma);

        if (counter >= line_count) {
          printf("Errore: counter %d supera line_count %d\n", counter, line_count); 
        }
        else if (strcmp(potential_evolution[counter], "0") == 0){
          if (counter == 0) {
            potential_zero = psi->psi->data[addr_rank0(psi->nsites, index)];
          }
          counter++; 
        }
        else if (strcmp(potential_evolution[counter], "0") != 0){
          if (contatore_update_potenziale == 1){
            if (strcmp(potential_evolution[counter], "t_0+") == 0){
              //psi->psi->data[addr_rank0(psi->nsites, index)] = potential_zero;
              psi->psi->data[addr_rank0(psi->nsites, index)] -= POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_1+") == 0){
              //psi->psi->data[addr_rank0(psi->nsites, index)] = potential_zero + POTENTIAL;
              psi->psi->data[addr_rank0(psi->nsites, index)] += POTENTIAL;
              counter++; 
            }
            else if (strcmp(potential_evolution[counter], "t_0-") == 0){
              //psi->psi->data[addr_rank0(psi->nsites, index)] = potential_zero;
              psi->psi->data[addr_rank0(psi->nsites, index)] -= -POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_1-") == 0){
              //psi->psi->data[addr_rank0(psi->nsites, index)] = potential_zero + POTENTIAL;
              psi->psi->data[addr_rank0(psi->nsites, index)] += -POTENTIAL;
              counter++; 
            }
            else {
              counter++;
            }
          }
          else if (contatore_update_potenziale == 2){
            if (strcmp(potential_evolution[counter], "t_1+") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_2+") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] += POTENTIAL;
              counter++; 
            }
            else if (strcmp(potential_evolution[counter], "t_1-") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= -POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_2-") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] += -POTENTIAL;
              counter++; 
            }
            else {
              counter++;
            }
           }
          else if (contatore_update_potenziale == 3){
            if (strcmp(potential_evolution[counter], "t_2+") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_3+") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] += POTENTIAL;
              counter++; 
            }
            else if (strcmp(potential_evolution[counter], "t_2-") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= -POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_3-") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] += -POTENTIAL;
              counter++; 
            }
            else {
              counter++;
            }
          }
          else if (contatore_update_potenziale == 4){
            if (strcmp(potential_evolution[counter], "t_3+") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= POTENTIAL;
              counter++;  
            } 
            else if (strcmp(potential_evolution[counter], "t_3-") == 0){
              psi->psi->data[addr_rank0(psi->nsites, index)] -= -POTENTIAL;
              counter++;  
            } 
            else {
              counter++;
            }
          }
        }
        
        
        else {
          printf("ERROR READING potential_evolution.txt \n");
	        exit(0);
        }

        //printf("row: %d psi: %f \n", counter, psi->psi->data[addr_rank0(psi->nsites, index)]);
        // LO COMMENTO PERCHE' NON SO SE SERVA
        /*
	      if (sigma) {
	        if (sigma > 0) {
	          psi_rho_set(psi, index, 0, sigma);
	          psi_rho_set(psi, index, 1, 0);
	       }
	       if (sigma < 0) {
	          psi_rho_set(psi, index, 0, 0);
	          psi_rho_set(psi, index, 1, -sigma);
	        }
	      }
        */ 
      }
      counter = counter + shift_z;
    }
    counter = counter + shift_y;
  }

  map_halo(map);
  

  return 0;
}

