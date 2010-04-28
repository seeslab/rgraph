#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include <gsl/gsl_rng.h>

#include "graph.h"
#include "tools.h"
#include "modules.h"
#include "layout.h"

/*
  Create a pointer to a node inside a simulation box
*/
struct node_box *
CreateNodeBox(struct node_gra *node)
{
  struct node_box *box;

  box = (struct node_box *)malloc(sizeof(struct node_box));
  box->ref = node;
  box->next = NULL;

  return box;
}

/*
  Add a node to a simulation box
*/
void AddNodeToBox(struct node_gra *node, struct node_box *box)
{
  struct node_box *temp = box->next;

  box->next = CreateNodeBox(node);
  box->next->next = temp;
}

/*
  Remove a node from a simulation box
*/
void RemoveNodeFromBox(struct node_gra *node, struct node_box *box)
{
  struct node_box *temp, *p = box;

  /* Find the node */
  while (p->next->ref != node)
    p = p->next;

  /* Remove the node */
  temp = p->next;
  p->next = temp->next;

  free(temp);
}

/*
  Free the memory occupied by a simulation box and its pointers
*/
void FreeBox(struct node_box *box)
{
  while (box->next != NULL)
    RemoveNodeFromBox(box->next->ref, box);
  free(box);
  return;
}

/*
  Print the coordinates of a node
*/
void
PrintNodeCoordinates(FILE *outFile, struct node_gra *net)
{
  while ((net = net->next) != NULL)
    fprintf(outFile, "%s %g %g\n", net->label, net->coorX, net->coorY);

  return;
}


/*
  Calculate the distance between circular areas for the component
  arrangement---returns -1 if target overlaps one or more other
  rectangles
*/
double CircleDistance(int target, double xtot, double ytot,
		      double *xce, double *yce, double *rad, int N)
{
  int i;
  double dis = 0.0;
  double overlap;

  /* Distance to the center */
  dis = sqrt((xce[target] - xtot) * (xce[target] - xtot) + 
	     (yce[target] - ytot) * (yce[target] - ytot));
  
  /* Overlap with other circles */
  for (i=0; i<N; i++) {
    if (i != target) {
      overlap = rad[target] + rad[i] -
	sqrt((xce[target] - xce[i]) * (xce[target] - xce[i]) + 
	     (yce[target] - yce[i]) * (yce[target] - yce[i]));
      if (overlap > 0.) {
	dis += 1e5 * overlap;
      }
    }
  }
  
  return dis;
}

/*
  Arrange isolated components
*/
void ArrangeComponents(struct node_gra *net, gsl_rng *gen)
{
  struct node_gra *cnet = NULL;
  struct node_gra *p = NULL, *cp = NULL;
  struct node_lis *node = NULL;
  struct group *part = NULL;
  struct group **glist = NULL;
  struct group *g = NULL;
  int nnod = CountNodes(net), ngroup;
  double *xce, *yce, *rad, *dxtot, *dytot;
  int i, target;
  double dis, totdis = 0.0, width, disi, disf, dx, dy;
  int test, maxsize = 0;
  double xtot, ytot;

  // Only do something if there is more than one component
  if ( IsGraphConnected(net) == 0 ) {
    cnet = CopyNetwork(net); // Create a copy of the network so the
			     // inGroup fields of the nodes are not
			     // overwritten
    part = ClustersPartition(cnet);
    ngroup = NNonEmptyGroups(part);
    
    // Allocate memory
    rad = allocate_d_vec(ngroup);
    xce = allocate_d_vec(ngroup);
    yce = allocate_d_vec(ngroup);
    dxtot = allocate_d_vec(ngroup);
    dytot = allocate_d_vec(ngroup);
    glist = (struct group **)malloc(ngroup * sizeof(struct group *));

    // Establish the "component circle"
    g = part;
    while (g->next != NULL) {
      g = g->next;
      glist[g->label] = g;
      dxtot[g->label] = dytot[g->label] = 0.0;

      // Find the center of the component and of the network
      node = g->nodeList;
      while(node->next != NULL) {
	node = node->next;
	xce[g->label] += node->ref->coorX;
	yce[g->label] += node->ref->coorY;
      }
      xce[g->label] /= (double)g->size;
      yce[g->label] /= (double)g->size;
      
      // Position of the center of the network
      if (g->size > maxsize) {
	maxsize = g->size;
	xtot = xce[g->label];
	ytot = yce[g->label];
      }

      // Find the radius of the circle
      rad[g->label] = 0.0;
      node = g->nodeList;
      while(node->next != NULL) {
	node = node->next;
	dis = sqrt((xce[g->label] - node->ref->coorX) *
		   (xce[g->label] - node->ref->coorX) +
		   (yce[g->label] - node->ref->coorY) *
		   (yce[g->label] - node->ref->coorY));
	if (dis > rad[g->label]) {
	  rad[g->label] = dis;
	}
      }
      rad[g->label] *= 1.02;
    }

    // Move the circles around
    width = 1.;
    while (width > 1.0e-2) {

      for (i=0; i<100*ngroup; i++) {
	// Select a group
	target = floor((double)ngroup * gsl_rng_uniform(gen));

	// Calculate the initial distance
	disi = CircleDistance(target, xtot, ytot,
			      xce, yce, rad, ngroup);
      
	// Move the group
	dx = width * 2. * (gsl_rng_uniform(gen) - .5);
	dy = width * 2. * (gsl_rng_uniform(gen) - .5);

	xce[target] += dx;
	yce[target] += dy;
      
	// Calculate the final distance
	disf = CircleDistance(target, xtot, ytot,
			      xce, yce, rad, ngroup);

	// Undo the movement if df>di
	if ( disf > disi ) {
	  xce[target] -= dx;
	  yce[target] -= dy;
	}
	else { // update total displacement
	  dxtot[target] += dx;
	  dytot[target] += dy;
	  totdis += disf - disi;
	}
      }
      // Decrease the amplitude of the movements;
      width *= 0.9;
    }

    // Reassing the coordinates of the nodes
    p = net;
    cp = cnet;
    while (p->next != NULL) {
      p = p->next;
      cp = cp->next;
      p->coorX += dxtot[cp->inGroup];
      p->coorY += dytot[cp->inGroup];
    }
	

    // Free memory
    free(xce);
    free(yce);
    free(rad);
    free(glist);
    RemoveGraph(cnet);
    RemovePartition(part);
  }
}

// Sets the X and Y coordinates of all the nodes between 0 and 1
void NormalizeCoordinates(struct node_gra *net)
{
  int nnod = CountNodes(net);
  double xmin = 1e100, xmax = -1e100;
  double ymin = 1e100, ymax = -1e100;
  struct node_gra *temp = NULL;

  // Find the minimum and maximum values of X and Y
  temp = net;
  while (temp->next != NULL) {
    temp = temp->next;

    if (temp->coorX < xmin)
      xmin = temp->coorX;
    if (temp->coorX > xmax)
      xmax = temp->coorX;
    if (temp->coorY < ymin)
      ymin = temp->coorY;
    if (temp->coorY > ymax)
      ymax = temp->coorY;
  }

  // Normalize the coordinates
  temp = net;
  while (temp->next != NULL) {
    temp = temp->next;

    temp->coorX = (temp->coorX - xmin) / (xmax - xmin);
    temp->coorY = (temp->coorY - ymin) / (ymax - ymin);
  }
}


// Sets the X and Y coordinates of all the nodes between 0 and 1
void NormalizeCoordinates3D(struct node_gra *net)
{
  int nnod = CountNodes(net);
  double xmin = 1e100, xmax = -1e100;
  double ymin = 1e100, ymax = -1e100;
  double zmin = 1e100, zmax = -1e100;
  struct node_gra *temp = NULL;

  // Find the minimum and maximum values of X and Y
  temp = net;
  while (temp->next != NULL) {
    temp = temp->next;

    if (temp->coorX < xmin)
      xmin = temp->coorX;
    if (temp->coorX > xmax)
      xmax = temp->coorX;
    if (temp->coorY < ymin)
      ymin = temp->coorY;
    if (temp->coorY > ymax)
      ymax = temp->coorY;
    if (temp->coorZ < zmin)
      zmin = temp->coorZ;
    if (temp->coorZ > zmax)
      zmax = temp->coorZ;
  }

  // Normalize the coordinates
  temp = net;
  while (temp->next != NULL) {
    temp = temp->next;

    temp->coorX = (temp->coorX - xmin) / (xmax - xmin);
    temp->coorY = (temp->coorY - ymin) / (ymax - ymin);
    temp->coorZ = (temp->coorZ - zmin) / (zmax - zmin);
  }
}


// Given the coordinates of a node, and the size of the system, find
// out in which box a node is
int AssignBox(struct node_gra *node, int nbox,
	      double xmin, double xmax, double ymin, double ymax)
{
  double normx, normy;
  int intx, inty;

  // Make the limits slightly larger than they are
  xmin = xmin - (xmax - xmin) * 1.e-5;
  xmax = xmax + (xmax - xmin) * 1.e-5;
  ymin = ymin - (ymax - ymin) * 1.e-5;
  ymax = ymax + (ymax - ymin) * 1.e-5;

  // Determine row and column
  normx = (node->coorX - xmin) / (xmax - xmin);
  normy = (node->coorY - ymin) / (ymax - ymin);
  intx = floor(normx * (double)nbox);
  inty = floor(normy * (double)nbox);

  return inty * nbox + intx;
}


// Calculates the harmonic and electrostatic potential energy of the
// network, using only the neighboring nodes for the calculation of
// the electrostatic interaction.
//
// CAUTION: Only works for undirected networks
double PotentialEnergy(struct node_gra *net,
		       struct node_box *box_list[], int boxprow)
{
  double E = 0.0;
  struct node_lis *nei;
  struct node_gra *n1 = NULL, *n2 = NULL;
  double hookek;
  double dis;
  double dx, dy;
  double Ke = 1. / (double)CountNodes(net);
  int nnod = CountNodes(net);
  int box, nei_box, nbox;
  struct node_box *nb1= NULL, *nb2= NULL;
  int row, col;
  int condition;
  double epsilon = sqrt(nnod) / 1.0e4;
  double xmin, xmax, ymin, ymax, cutoff;

  // Calculate the cutoff
  xmin = ymin = 1.e5;
  xmax = ymax = -1.e5;
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;
    if (n1->coorX < xmin)
      xmin = n1->coorX;
    if (n1->coorY < ymin)
      ymin = n1->coorY;
    if (n1->coorX > xmax)
      xmax = n1->coorX;
    if (n1->coorY > ymax)
      ymax = n1->coorY;
  }
/*   cutoff = min(xmax-xmin, ymax-ymin); */
  cutoff = (xmax - xmin + ymax - ymin) / 2.;
  cutoff /= (double)boxprow;

  // Harmonic interaction with neighbors
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;

    nei = n1->neig;
    while (nei->next != NULL){
      nei = nei->next;
      n2 = nei->ref;

      if (n1->num > n2->num) {
	hookek = nei->weight;
	dx = n2->coorX - n1->coorX;
	dy = n2->coorY - n1->coorY;
	dis = sqrt( dx * dx + dy * dy );

	E += hookek * dis * dis / 2.0;
      }
    }
  }

  // Electrostatic interaction with all other nodes in neighboring
  // boxes
  nbox = boxprow * boxprow;
  for (box=0; box<nbox; box++) {
    
    // Neighboring boxes
    for (row=0; row<=1; row++) {
      for (col=-1; col<=1; col++) {
	nei_box = box + row * boxprow + col;
	  
	// Do not consider some "neighboring" boxes that are, in
	// reality, not neighbors of the current box, and avoid double
	// counting of pairs of boxes.
	condition = (nei_box < box) || // box to the left
	  (col == 1 && nei_box % boxprow == 0) || // rightmost col
	  (col == -1 && box % boxprow == 0) || // leftmost col
	  (row == 1 && box >= nbox-boxprow); // bottom row

	// Calculate forces
	if ( !condition ) {

	  nb1 = box_list[box];
	  while (nb1->next != NULL) {
	    nb1 = nb1->next;
	    n1 = nb1->ref;

	    nb2 = box_list[nei_box];
	    while (nb2->next != NULL) {
	      nb2 = nb2->next;
	      n2 = nb2->ref;
	  
	      if ( (box != nei_box) ||
		   ( (box == nei_box) && (n1->num > n2->num) ) ) {

		dx = n2->coorX - n1->coorX;
		dy = n2->coorY - n1->coorY;
		dis = sqrt( dx * dx + dy * dy );
		if (dis < cutoff || boxprow <= 2) {
		  E += Ke / sqrt(dis*dis + epsilon);
		}
	      }
	    }
	  } // Energy calculated for this pair of boxes
	}
      }
    }  // End of loop over neighboring boxes

  } // End of loop over boxes

  return E;
}


// Calculates the harmonic and electrostatic forces felt by a node,
// using only the neighboring nodes for the calculation of the
// electrostatic interaction.
//
// CAUTION: Only works for undirected networks
void
CalculateNodeForces(struct node_gra *net, double Fx[], double Fy[],
		    struct node_box *box_list[], int boxprow)
{
  struct node_lis *nei;
  struct node_gra *n1 = NULL, *n2 = NULL;
  double hookek;
  double dis;
  double dx, dy, Fmod;
  double Ke = 1. / (double)CountNodes(net);
  int nnod = CountNodes(net);
  int box, nei_box, nbox;
  struct node_box *nb1= NULL, *nb2= NULL;
  int row, col;
  int condition;
  double epsilon = sqrt(nnod) / 1.0e4;
  double xmin, xmax, ymin, ymax, cutoff;

  // Set all forces to 0.0 and calculate the size of the boxes
  xmin = ymin = 1.e5;
  xmax = ymax = -1.e5;
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;
    Fx[n1->trans] = Fy[n1->trans] = 0.0;
    if (n1->coorX < xmin)
      xmin = n1->coorX;
    if (n1->coorY < ymin)
      ymin = n1->coorY;
    if (n1->coorX > xmax)
      xmax = n1->coorX;
    if (n1->coorY > ymax)
      ymax = n1->coorY;
  }
/*   cutoff = min(xmax-xmin, ymax-ymin); */
  cutoff = (xmax - xmin + ymax - ymin) / 2.;
  cutoff /= (double)boxprow;

  // Harmonic interaction with neighbors
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;

    nei = n1->neig;
    while (nei->next != NULL){
      nei = nei->next;
      n2 = nei->ref;

      if (n1->num > n2->num) {
	hookek = nei->weight;
	dx = n2->coorX - n1->coorX;
	dy = n2->coorY - n1->coorY;
	dis = sqrt( dx * dx + dy * dy );
	Fmod = hookek * dis;
	
	Fx[n1->trans] += Fmod * dx / dis;
	Fy[n1->trans] += Fmod * dy / dis;
	Fx[n2->trans] -= Fmod * dx / dis;
	Fy[n2->trans] -= Fmod * dy / dis;
      }
    }
  }

  // Electrostatic interaction with all other nodes in neighboring
  // boxes
  nbox = boxprow * boxprow;
  for (box=0; box<nbox; box++) {
    
    // Neighboring boxes
    for (row=0; row<=1; row++) {
      for (col=-1; col<=1; col++) {
	nei_box = box + row * boxprow + col;
	  
	// Do not consider some "neighboring" boxes that are, in
	// reality, not neighbors of the current box, and avoid double
	// counting of pairs of boxes.
	condition = (nei_box < box) || // box to the left
	  (col == 1 && nei_box % boxprow == 0) || // rightmost col
	  (col == -1 && box % boxprow == 0) || // leftmost col
	  (row == 1 && box >= nbox-boxprow); // bottom row

	// Calculate forces
	if ( !condition ) {

	  nb1 = box_list[box];
	  while (nb1->next != NULL) {
	    nb1 = nb1->next;
	    n1 = nb1->ref;

	    nb2 = box_list[nei_box];
	    while (nb2->next != NULL) {
	      nb2 = nb2->next;
	      n2 = nb2->ref;
	  
	      if ( (box != nei_box) ||
		   ( (box == nei_box) && (n1->num > n2->num) ) ) {

		dx = n2->coorX - n1->coorX;
		dy = n2->coorY - n1->coorY;
		dis = sqrt( dx * dx + dy * dy );
		if (dis < cutoff || boxprow <= 2) {
		  Fmod = -Ke / (dis*dis + epsilon);

		  Fx[n1->trans] += Fmod * dx / pow(dis, 1.);
		  Fy[n1->trans] += Fmod * dy / pow(dis, 1.);
		  Fx[n2->trans] -= Fmod * dx / pow(dis, 1.);
		  Fy[n2->trans] -= Fmod * dy / pow(dis, 1.);
		}
	      }
	    }
	  } // Forces calculated for this pair of boxes
	}
      }
    }  // End of loop over neighboring boxes

  } // End of loop over boxes

  return;
}

// Layout the network using the spring+charge+drag forces and the
// Verlet algorithm
void
MDGraphLayout(struct node_gra *net, double drag, double dt,
	      int nsteps, gsl_rng *gen, int nbox)
{
  double *Fx = NULL, *Fy = NULL;
  double *prev_xco = NULL, *prev_yco = NULL;
  double temp_x, temp_y;
  int nnod = CountNodes(net);
  int ncoun = 0;
  struct node_gra *p;
  int condition;
  int i;
  double vx, vy, T=0.; // Kinetic energy
  double U; // potential energy
  int nb2;
  struct node_box **box_list;
  int *oldbox, newbox;
  double xmin, xmax, ymin, ymax;

  // Determine the number of boxes
  if (nbox <= 0) {
    nbox = floor(sqrt((double)(nnod*nnod) / 
		      (10.*(double)TotalNLinks(net, 1))));
    if (nbox < 1)
      nbox = 1;
  }
  nb2 = nbox * nbox;
  fprintf(stderr, "# Using %dx%d boxes\n", nbox, nbox);

  // Allocate memory
  Fx = allocate_d_vec(nnod);
  Fy = allocate_d_vec(nnod);
  prev_xco = allocate_d_vec(nnod);
  prev_yco = allocate_d_vec(nnod);
  oldbox = allocate_i_vec(nnod);

  box_list = (struct node_box **)malloc(nb2*sizeof(struct node_box *));
  for (i=0; i<nb2; i++) {
    box_list[i] = CreateNodeBox(NULL);
  }
  
  // Set initial conditions and the translation for each node
  p = net;
  ncoun = 0;
  while (p->next != NULL) {
    p = p->next;
    p->coorX = gsl_rng_uniform(gen);
    p->coorY = gsl_rng_uniform(gen);
    p->trans = ncoun; // trans stores the translation of the node
    prev_xco[ncoun] = p->coorX; // Zero X velocity
    prev_yco[ncoun] = p->coorY; // Zero Y velocity
    oldbox[ncoun] = AssignBox(p, nbox, 0.0, 1.0, 0.0, 1.0);
    AddNodeToBox(p, box_list[oldbox[ncoun]]);
    ncoun++;
  }

  // Integrate the equations of motion
  for (i=0; i<nsteps; i++) {

    // Reset the size of the box
    xmin = ymin = 1.e50;
    xmax = ymax = -1.e50;
    
/*     // Calculate the potential energy at time t */
/*     U = PotentialEnergy(net, box_list, nbox); */

    // Calculate the forces
    CalculateNodeForces(net, Fx, Fy, box_list, nbox);
  
    T = 0.;
    // Update node's positions using Verlet and adding a drag term
    // -drag*v(t) to the force
    p = net;
    while (p->next != NULL) {
      p = p->next;
      
      temp_x = p->coorX;
      temp_y = p->coorY;
      
      p->coorX = (2. * temp_x -
		  prev_xco[p->trans] * (1. - drag * dt * .5) +
		  Fx[p->trans] * dt * dt) /
	(1. + drag * dt * .5);
      p->coorY = (2. * temp_y -
		  prev_yco[p->trans] * (1. - drag * dt * .5) +
		  Fy[p->trans] * dt * dt) /
	(1. + drag * dt * .5);

      // Re-check the size of the box
      if (p->coorX < xmin)
	xmin = p->coorX;
      if (p->coorY < ymin)
	ymin = p->coorY;
      if (p->coorX > xmax)
	xmax = p->coorX;
      if (p->coorY > ymax)
	ymax = p->coorY;

      // Calculate the kinetic energy at time t
      vx = (p->coorX - prev_xco[p->trans]) / (2. * dt);
      vy = (p->coorY - prev_yco[p->trans]) / (2. * dt);
      T += (vx*vx + vy*vy) / 2.;

      prev_xco[p->trans] = temp_x;
      prev_yco[p->trans] = temp_y;

    } // All positions have been updated


    // Move the nodes to new boxes if necessary
    p = net;
    while (p->next != NULL) {
      p = p->next;
      newbox = AssignBox(p, nbox, xmin, xmax, ymin, ymax);
      if (newbox != oldbox[p->trans]) {
	RemoveNodeFromBox(p, box_list[oldbox[p->trans]]);
	AddNodeToBox(p, box_list[newbox]);
	oldbox[p->trans] = newbox;
      }
    }

    // Screen output
    fprintf(stderr, "%d %g\n", i, T/(double)nnod);
/*     printf("%d %g %g %g\n", i, T/(double)nnod, */
/* 	   U/(double)nnod, */
/* 	   (T+U)/(double)nnod); */
/*     printf("%d\n", i); */

    // Check if nodes are still moving
    if ( T/(double)nnod < 1.e-5 ) {
      i = nsteps;                   // Stop the simulation
    }

  } // End of MD

  // Normalize the coordinates
  NormalizeCoordinates(net);
  ArrangeComponents(net, gen);
  NormalizeCoordinates(net);

  // Free memory
  free(Fx);
  free(Fy);
  free(prev_xco);
  free(prev_yco);
  free(oldbox);
  for (i=0; i<nb2; i++)
    FreeBox(box_list[i]);
  free(box_list);
}

// Given the coordinates of a node, and the size of the system, find
// out in which box a node is
int AssignBox3D(struct node_gra *node, int nbox,
		double xmin, double xmax,
		double ymin, double ymax,
		double zmin, double zmax)
{
  double normx, normy, normz;
  int intx, inty, intz;

  // Make the limits slightly larger than they are
  xmin = xmin - (xmax - xmin) * 1.e-5;
  xmax = xmax + (xmax - xmin) * 1.e-5;
  ymin = ymin - (ymax - ymin) * 1.e-5;
  ymax = ymax + (ymax - ymin) * 1.e-5;
  zmin = zmin - (zmax - zmin) * 1.e-5;
  zmax = zmax + (zmax - zmin) * 1.e-5;

  // Determine row and column
  normx = (node->coorX - xmin) / (xmax - xmin);
  normy = (node->coorY - ymin) / (ymax - ymin);
  normz = (node->coorZ - zmin) / (zmax - zmin);
  intx = floor(normx * (double)nbox);
  inty = floor(normy * (double)nbox);
  intz = floor(normz * (double)nbox);

  return intz * nbox * nbox + inty * nbox + intx;
}


// Calculates the harmonic and electrostatic potential energy of the
// network, using only the neighboring nodes for the calculation of
// the electrostatic interaction.
//
// CAUTION: Only works for undirected networks
double PotentialEnergy3D(struct node_gra *net,
			 struct node_box *box_list[], int boxprow)
{
  double E = 0.0;
  struct node_lis *nei;
  struct node_gra *n1 = NULL, *n2 = NULL;
  double hookek;
  double dis2;
  double dx, dy, dz;
  double Ke = 1. / (double)CountNodes(net);
  int nnod = CountNodes(net);
  int box, nei_box, nbox;
  struct node_box *nb1= NULL, *nb2= NULL;
  int row, col, pla;
  int condition;
  double epsilon = sqrt(nnod) / 1.0e4;
  double xmin, xmax, ymin, ymax, zmin, zmax, cutoff, cutoff2;

  // Calculate the cutoff
  xmin = ymin = zmin = 1.e5;
  xmax = ymax = zmax = -1.e5;
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;
    if (n1->coorX < xmin)
      xmin = n1->coorX;
    if (n1->coorY < ymin)
      ymin = n1->coorY;
    if (n1->coorZ < zmin)
      zmin = n1->coorZ;
    if (n1->coorX > xmax)
      xmax = n1->coorX;
    if (n1->coorY > ymax)
      ymax = n1->coorY;
    if (n1->coorZ > zmax)
      zmax = n1->coorZ;
  }
  cutoff = (xmax - xmin + ymax - ymin + zmax - zmin) / 3.;
  cutoff /= (double)boxprow;
  cutoff2 = cutoff * cutoff;

  // Harmonic interaction with neighbors
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;

    nei = n1->neig;
    while (nei->next != NULL){
      nei = nei->next;
      n2 = nei->ref;

      if (n1->num > n2->num) {
	hookek = nei->weight;
	dx = n2->coorX - n1->coorX;
	dy = n2->coorY - n1->coorY;
	dz = n2->coorZ - n1->coorZ;
	dis2 = dx * dx + dy * dy + dz * dz;

	E += hookek * dis2 / 2.0;
      }
    }
  }

  // Electrostatic interaction with all other nodes in neighboring
  // boxes
  nbox = boxprow * boxprow;
  for (box=0; box<nbox; box++) {
    
    // Neighboring boxes
    for (row=-1; row<=1; row++) {
      for (col=-1; col<=1; col++) {
	for (pla=-1; pla<=1; pla++) {
	  nei_box = box +
	    pla * boxprow * boxprow + row * boxprow + col;
	  
	  // Do not consider some "neighboring" boxes that are, in
	  // reality, not neighbors of the current box, and avoid
	  // double counting of pairs of boxes.
	  condition =
	    // avoid double counting
	    (nei_box < box) ||
	    // rightmost col
	    (col == 1 && nei_box % boxprow == 0) ||
	    // leftmost col
	    (col == -1 && box % boxprow == 0) ||
	    // bottom row
	    (row == 1 && nei_box % (boxprow*boxprow) < boxprow) ||
	    // top row
	    (row == -1 && box % (boxprow*boxprow) < boxprow) ||
	    // last plane
	    (pla == 1 && box >= nbox-boxprow*boxprow) ||
	    // first plane
	    (pla == -1 && box < boxprow*boxprow);

	  // Calculate forces
	  if ( !condition ) {

	    nb1 = box_list[box];
	    while (nb1->next != NULL) {
	      nb1 = nb1->next;
	      n1 = nb1->ref;
	      
	      nb2 = box_list[nei_box];
	      while (nb2->next != NULL) {
		nb2 = nb2->next;
		n2 = nb2->ref;
		
		if ( (box != nei_box) ||
		     ( (box == nei_box) && (n1->num > n2->num) ) ) {
		  
		  dx = n2->coorX - n1->coorX;
		  dy = n2->coorY - n1->coorY;
		  dz = n2->coorZ - n1->coorZ;
		  dis2 = dx * dx + dy * dy + dz * dz;
		  if (dis2 < cutoff2 || boxprow <= 2) {
		    E += Ke / sqrt(dis2 + epsilon);
		  }
		}
	      }
	    } // Energy calculated for this pair of boxes
	  }
	}
      }
    }  // End of loop over neighboring boxes

  } // End of loop over boxes

  return E;
}


// Calculates the harmonic and electrostatic forces felt by a node,
// using only the neighboring nodes for the calculation of the
// electrostatic interaction.
//
// CAUTION: Only works for undirected networks
void
CalculateNodeForces3D(struct node_gra *net,
		      double Fx[], double Fy[], double Fz[],
		      struct node_box *box_list[], int boxprow)
{
  struct node_lis *nei;
  struct node_gra *n1 = NULL, *n2 = NULL;
  double hookek;
  double dis;
  double dx, dy, dz, Fmod;
  double Ke = 1. / (double)CountNodes(net);
  int nnod = CountNodes(net);
  int box, nei_box, nbox;
  struct node_box *nb1= NULL, *nb2= NULL;
  int row, col, pla;
  int condition;
  double epsilon = sqrt(nnod) / 1.0e4;
  double xmin, xmax, ymin, ymax, zmin, zmax, cutoff;

  // Set all forces to 0.0 and calculate the size of the boxes
  xmin = ymin = zmin = 1.e5;
  xmax = ymax = zmax = -1.e5;
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;
    Fx[n1->trans] = Fy[n1->trans] = Fz[n1->trans] = 0.0;
    if (n1->coorX < xmin)
      xmin = n1->coorX;
    if (n1->coorY < ymin)
      ymin = n1->coorY;
    if (n1->coorZ < zmin)
      zmin = n1->coorZ;
    if (n1->coorX > xmax)
      xmax = n1->coorX;
    if (n1->coorY > ymax)
      ymax = n1->coorY;
    if (n1->coorZ > zmax)
      zmax = n1->coorZ;
  }
  cutoff = (xmax - xmin + ymax - ymin + zmax - zmin) / 3.;
  cutoff /= (double)boxprow;

  // Harmonic interaction with neighbors
  n1 = net;
  while (n1->next != NULL) {
    n1 = n1->next;

    nei = n1->neig;
    while (nei->next != NULL){
      nei = nei->next;
      n2 = nei->ref;

      if (n1->num > n2->num) {
	hookek = nei->weight;
	dx = n2->coorX - n1->coorX;
	dy = n2->coorY - n1->coorY;
	dz = n2->coorZ - n1->coorZ;
	dis = sqrt( dx * dx + dy * dy + dz * dz);
	Fmod = hookek * dis;
	
	Fx[n1->trans] += Fmod * dx / dis;
	Fy[n1->trans] += Fmod * dy / dis;
	Fz[n1->trans] += Fmod * dz / dis;
	Fx[n2->trans] -= Fmod * dx / dis;
	Fy[n2->trans] -= Fmod * dy / dis;
	Fz[n2->trans] -= Fmod * dz / dis;
      }
    }
  }

  // Electrostatic interaction with all other nodes in neighboring
  // boxes
  nbox = boxprow * boxprow * boxprow;
  for (box=0; box<nbox; box++) {
    
    // Neighboring boxes
    for (row=-1; row<=1; row++) {
      for (col=-1; col<=1; col++) {
	for (pla=-1; pla<=1; pla++) {
	  nei_box = box +
	    pla * boxprow * boxprow + row * boxprow + col;
	  
	  // Do not consider some "neighboring" boxes that are, in
	  // reality, not neighbors of the current box, and avoid
	  // double counting of pairs of boxes.
	  condition =
	    // avoid double counting
	    (nei_box < box) ||
	    // rightmost col
	    (col == 1 && nei_box % boxprow == 0) ||
	    // leftmost col
	    (col == -1 && box % boxprow == 0) ||
	    // bottom row
	    (row == 1 && nei_box % (boxprow*boxprow) < boxprow) ||
	    // top row
	    (row == -1 && box % (boxprow*boxprow) < boxprow) ||
	    // last plane
	    (pla == 1 && box >= nbox-boxprow*boxprow) ||
	    // first plane
	    (pla == -1 && box < boxprow*boxprow);

	  // Calculate forces
	  if ( !condition ) {

	    nb1 = box_list[box];
	    while (nb1->next != NULL) {
	      nb1 = nb1->next;
	      n1 = nb1->ref;

	      nb2 = box_list[nei_box];
	      while (nb2->next != NULL) {
		nb2 = nb2->next;
		n2 = nb2->ref;
	  
		if ( (box != nei_box) ||
		     ( (box == nei_box) && (n1->num > n2->num) ) ) {

		  dx = n2->coorX - n1->coorX;
		  dy = n2->coorY - n1->coorY;
		  dz = n2->coorZ - n1->coorZ;
		  dis = sqrt( dx * dx + dy * dy + dz * dz);
		  if (dis < cutoff || boxprow <= 2) {
		    Fmod = -Ke / (dis*dis + epsilon);

		    Fx[n1->trans] += Fmod * dx / dis;
		    Fy[n1->trans] += Fmod * dy / dis;
		    Fz[n1->trans] += Fmod * dz / dis;
		    Fx[n2->trans] -= Fmod * dx / dis;
		    Fy[n2->trans] -= Fmod * dy / dis;
		    Fz[n2->trans] -= Fmod * dz / dis;
		  }
		}
	      }
	    } // Forces calculated for this pair of boxes
	  }
	}
      }
    }  // End of loop over neighboring boxes

  } // End of loop over boxes

  return;
}

// Layout the network in 3D using the spring+charge+drag forces and
// the Verlet algorithm
void
MDGraphLayout3D(struct node_gra *net, double drag, double dt,
		int nsteps, gsl_rng *gen, int nbox)
{
  double *Fx = NULL, *Fy = NULL, *Fz = NULL;
  double *prev_xco = NULL, *prev_yco = NULL, *prev_zco = NULL;
  double temp_x, temp_y, temp_z;
  int nnod = CountNodes(net);
  int ncoun = 0;
  struct node_gra *p;
  int condition;
  int i;
  double vx, vy, vz, T=0.; // Kinetic energy
  double U; // potential energy
  int nb3;
  struct node_box **box_list;
  int *oldbox, newbox;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  // Determine the number of boxes
  if (nbox <= 0) {
    nbox = floor(pow((double)(nnod*nnod) / 
		     (10.*(double)TotalNLinks(net, 1)), 1./3.));
    if (nbox < 1)
      nbox = 1;
  }
  nb3 = nbox * nbox * nbox;
  printf("# Using %dx%dx%d boxes\n", nbox, nbox, nbox);

  // Allocate memory
  Fx = allocate_d_vec(nnod);
  Fy = allocate_d_vec(nnod);
  Fz = allocate_d_vec(nnod);
  prev_xco = allocate_d_vec(nnod);
  prev_yco = allocate_d_vec(nnod);
  prev_zco = allocate_d_vec(nnod);
  oldbox = allocate_i_vec(nnod);

  box_list = (struct node_box **)malloc(nb3*sizeof(struct node_box *));
  for (i=0; i<nb3; i++) {
    box_list[i] = CreateNodeBox(NULL);
  }
  
  // Set initial conditions and the translation for each node
  p = net;
  ncoun = 0;
  while (p->next != NULL) {
    p = p->next;
    p->coorX = gsl_rng_uniform(gen);
    p->coorY = gsl_rng_uniform(gen);
    p->coorZ = gsl_rng_uniform(gen);
    p->trans = ncoun; // trans stores the translation of the node
    prev_xco[ncoun] = p->coorX; // Zero X velocity
    prev_yco[ncoun] = p->coorY; // Zero Y velocity
    prev_zco[ncoun] = p->coorZ; // Zero Z velocity
    oldbox[ncoun] = AssignBox3D(p, nbox, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    AddNodeToBox(p, box_list[oldbox[ncoun]]);
    ncoun++;
  }

  // Integrate the equations of motion
  for (i=0; i<nsteps; i++) {

    // Reset the size of the box
    xmin = ymin = zmin = 1.e50;
    xmax = ymax = zmax = -1.e50;
    
/*     // Calculate the potential energy at time t */
/*     U = PotentialEnergy3D(net, box_list, nbox); */

    // Calculate the forces
    CalculateNodeForces3D(net, Fx, Fy, Fz, box_list, nbox);
  
    T = 0.;
    // Update node's positions using Verlet and adding a drag term
    // -drag*v(t) to the force
    p = net;
    while (p->next != NULL) {
      p = p->next;
      
      temp_x = p->coorX;
      temp_y = p->coorY;
      temp_z = p->coorZ;
      
      p->coorX = (2. * temp_x -
		  prev_xco[p->trans] * (1. - drag * dt * .5) +
		  Fx[p->trans] * dt * dt) /
	(1. + drag * dt * .5);
      p->coorY = (2. * temp_y -
		  prev_yco[p->trans] * (1. - drag * dt * .5) +
		  Fy[p->trans] * dt * dt) /
	(1. + drag * dt * .5);
      p->coorZ = (2. * temp_z -
		  prev_zco[p->trans] * (1. - drag * dt * .5) +
		  Fz[p->trans] * dt * dt) /
	(1. + drag * dt * .5);

      // Re-check the size of the box
      if (p->coorX < xmin)
	xmin = p->coorX;
      if (p->coorY < ymin)
	ymin = p->coorY;
      if (p->coorZ < zmin)
	zmin = p->coorZ;
      if (p->coorX > xmax)
	xmax = p->coorX;
      if (p->coorY > ymax)
	ymax = p->coorY;
      if (p->coorZ > zmax)
	zmax = p->coorZ;

      // Calculate the kinetic energy at time t
      vx = (p->coorX - prev_xco[p->trans]) / (2. * dt);
      vy = (p->coorY - prev_yco[p->trans]) / (2. * dt);
      vz = (p->coorZ - prev_zco[p->trans]) / (2. * dt);
      T += (vx*vx + vy*vy + vz*vz) / 2.;

      prev_xco[p->trans] = temp_x;
      prev_yco[p->trans] = temp_y;
      prev_zco[p->trans] = temp_z;

    } // All positions have been updated

    // Move the nodes to new boxes if necessary
    p = net;
    while (p->next != NULL) {
      p = p->next;
      newbox = AssignBox3D(p, nbox,
			   xmin, xmax, ymin, ymax, zmin, zmax);
      if (newbox != oldbox[p->trans]) {
	RemoveNodeFromBox(p, box_list[oldbox[p->trans]]);
	AddNodeToBox(p, box_list[newbox]);
	oldbox[p->trans] = newbox;
      }
    }

    // Screen output
    printf("%d %g\n", i, T/(double)nnod);
/*     printf("%d %g %g %g\n", i, T/(double)nnod, */
/* 	   U/(double)nnod, */
/* 	   (T+U)/(double)nnod); */
/*     printf("%d\n", i); */

    // Check if nodes are still moving
    if ( T/(double)nnod < 1.e-5 ) {
      i = nsteps;                   // Stop the simulation
    }

  } // End of MD

  // Normalize the coordinates
  NormalizeCoordinates3D(net);

  // Free memory
  free(Fx);
  free(Fy);
  free(Fz);
  free(prev_xco);
  free(prev_yco);
  free(prev_zco);
  free(oldbox);
  for (i=0; i<nb3; i++)
    FreeBox(box_list[i]);
  free(box_list);
}


// Layout the network using the spring+charge+drag forces and the
// Verlet algorithm. Use 3D first, and progressively pull the nodes
// towards the z=0.5 plane.
void
MDGraphLayout2Dp(struct node_gra *net, double drag, double dt,
		 int nsteps, gsl_rng *gen, int nbox)
{
  double *Fx = NULL, *Fy = NULL, *Fz = NULL;
  double *prev_xco = NULL, *prev_yco = NULL, *prev_zco = NULL;
  double temp_x, temp_y, temp_z;
  int nnod = CountNodes(net);
  int ncoun = 0;
  struct node_gra *p;
  int condition;
  int i;
  double vx, vy, vz, T=0.; // Kinetic energy
  double U; // potential energy
  int nb3;
  struct node_box **box_list;
  int *oldbox, newbox;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double kz = 0.0;
  int kz_sw = 0;

  // Determine the number of boxes
  if (nbox <= 0) {
    nbox = floor(pow((double)(nnod*nnod) / 
		     (10.*(double)TotalNLinks(net, 1)), 1./3.));
    if (nbox < 1)
      nbox = 1;
  }
  nb3 = nbox * nbox * nbox;
  printf("# Using %dx%dx%d boxes\n", nbox, nbox, nbox);

  // Allocate memory
  Fx = allocate_d_vec(nnod);
  Fy = allocate_d_vec(nnod);
  Fz = allocate_d_vec(nnod);
  prev_xco = allocate_d_vec(nnod);
  prev_yco = allocate_d_vec(nnod);
  prev_zco = allocate_d_vec(nnod);
  oldbox = allocate_i_vec(nnod);

  box_list = (struct node_box **)malloc(nb3*sizeof(struct node_box *));
  for (i=0; i<nb3; i++) {
    box_list[i] = CreateNodeBox(NULL);
  }
  
  // Set initial conditions and the translation for each node
  p = net;
  ncoun = 0;
  while (p->next != NULL) {
    p = p->next;
    p->coorX = gsl_rng_uniform(gen);
    p->coorY = gsl_rng_uniform(gen);
    p->coorZ = gsl_rng_uniform(gen);
    p->trans = ncoun; // trans stores the translation of the node
    prev_xco[ncoun] = p->coorX; // Zero X velocity
    prev_yco[ncoun] = p->coorY; // Zero Y velocity
    prev_zco[ncoun] = p->coorZ; // Zero Z velocity
    oldbox[ncoun] = AssignBox3D(p, nbox, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    AddNodeToBox(p, box_list[oldbox[ncoun]]);
    ncoun++;
  }

  // Integrate the equations of motion
  for (i=0; i<nsteps; i++) {

    // Reset the size of the box
    xmin = ymin = zmin = 1.e50;
    xmax = ymax = zmax = -1.e50;
    
/*     // Calculate the potential energy at time t */
/*     U = PotentialEnergy3D(net, box_list, nbox); */

    // Calculate the forces
    CalculateNodeForces3D(net, Fx, Fy, Fz, box_list, nbox);
  
    T = 0.;
    // Update node's positions using Verlet and adding a drag term
    // -drag*v(t) to the force
    p = net;
    while (p->next != NULL) {
      p = p->next;
      
      temp_x = p->coorX;
      temp_y = p->coorY;
      temp_z = p->coorZ;
      
      p->coorX = (2. * temp_x -
		  prev_xco[p->trans] * (1. - drag * dt * .5) +
		  Fx[p->trans] * dt * dt) /
	(1. + drag * dt * .5);
      p->coorY = (2. * temp_y -
		  prev_yco[p->trans] * (1. - drag * dt * .5) +
		  Fy[p->trans] * dt * dt) /
	(1. + drag * dt * .5);
      p->coorZ = (2. * temp_z -
		  prev_zco[p->trans] * (1. - drag * dt * .5) +
		  (Fz[p->trans] - kz * (temp_z - .5)) * dt * dt) /
	(1. + drag * dt * .5);

      // Re-check the size of the box
      if (p->coorX < xmin)
	xmin = p->coorX;
      if (p->coorY < ymin)
	ymin = p->coorY;
      if (p->coorZ < zmin)
	zmin = p->coorZ;
      if (p->coorX > xmax)
	xmax = p->coorX;
      if (p->coorY > ymax)
	ymax = p->coorY;
      if (p->coorZ > zmax)
	zmax = p->coorZ;

      // Calculate the kinetic energy at time t
      vx = (p->coorX - prev_xco[p->trans]) / (2. * dt);
      vy = (p->coorY - prev_yco[p->trans]) / (2. * dt);
      vz = (p->coorZ - prev_zco[p->trans]) / (2. * dt);
      T += (vx*vx + vy*vy + vz*vz) / 2.;

      prev_xco[p->trans] = temp_x;
      prev_yco[p->trans] = temp_y;
      prev_zco[p->trans] = temp_z;

    } // All positions have been updated

    // Move the nodes to new boxes if necessary
    p = net;
    while (p->next != NULL) {
      p = p->next;
      newbox = AssignBox3D(p, nbox,
			   xmin, xmax, ymin, ymax, zmin, zmax);
      if (newbox != oldbox[p->trans]) {
	RemoveNodeFromBox(p, box_list[oldbox[p->trans]]);
	AddNodeToBox(p, box_list[newbox]);
	oldbox[p->trans] = newbox;
      }
    }

    // Screen output
    printf("%d %g\n", i, T/(double)nnod);
/*     printf("%d %g %g %g\n", i, T/(double)nnod, */
/* 	   U/(double)nnod, */
/* 	   (T+U)/(double)nnod); */
/*     printf("%d\n", i); */

    // Increase the spring in the z direction
    if ( kz_sw > 0 ) {
      kz += dt;
    }

    // Check if nodes are still moving
    if ( (T/(double)nnod < 1.e-3) && (kz_sw == 0) ) {
      kz_sw = 1;         // Switch on the spring in the z direction
    }
    if ( T/(double)nnod < 1.e-5 ) {
      i = nsteps;                   // Stop the simulation
    }

  } // End of MD

  // Normalize the coordinates
  NormalizeCoordinates(net);
  ArrangeComponents(net, gen);
  NormalizeCoordinates(net);

  // Free memory
  free(Fx);
  free(Fy);
  free(Fz);
  free(prev_xco);
  free(prev_yco);
  free(prev_zco);
  free(oldbox);
  for (i=0; i<nb3; i++)
    FreeBox(box_list[i]);
  free(box_list);
}


// Assigns coordinates to the nodes in net_des from those in
// net_src. The program will return -1 and complain if a node from
// net_des is NOT found in net_src.
/* int CopyCoords(struct node_gra *net_src, struct node_gra *net_des) */
/* { */
/*   struct node_gra *p = NULL; */
/*   struct node_gra *nlist[maxim_int]; */

/*   // Map the src network into the nlist array for fast access */
/*   p = net_src; */
/*   while(p->next != NULL) { */
/*     p = p->next; */
/*     nlist[p->num] = p; */
/*   } */

/*   // Copy the coordinates */
/*   p = net_des; */
/*   while(p->next != NULL) { */
/*     p = p->next; */
/*     if (nlist[p->num] == NULL) { */
/*       printf("ERROR in CopyCoords! Quitting\n"); */
/*       return -1; */
/*     } */
/*     p->coorX = nlist[p->num]->coorX; */
/*     p->coorY = nlist[p->num]->coorY; */
/*   } */

/*   // Ended successfully */
/*   return 0; */
/* } */
