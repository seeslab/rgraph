#ifndef RGRAPH_VISUAL_H
#define RGRAPH_VISUAL_H 1

#include <gsl/gsl_rng.h>

/*
  -----------------------------------------------------------------------------
  Structure for a faster electrostatic force and energy calculation
  -----------------------------------------------------------------------------
*/
struct node_box{
  struct node_gra *ref;
  struct node_box *next;   /* next element in the list */
};


/*
  -----------------------------------------------------------------------------
  Box operations
  -----------------------------------------------------------------------------
*/
struct node_box *CreateNodeBox(struct node_gra *node);
void AddNodeToBox(struct node_gra *node, struct node_box *box);
void RemoveNodeFromBox(struct node_gra *node, struct node_box *box);
int AssignBox3D(struct node_gra *node, int nbox,
		double xmin, double xmax,
		double ymin, double ymax,
		double zmin, double zmax);
void FreeBox(struct node_box *box);

/*
  -----------------------------------------------------------------------------
  Distance, force, energy
  -----------------------------------------------------------------------------
*/
double CircleDistance(int target, double xtot, double ytot,
		      double *xce, double *yce, double *rad, int N);
void ArrangeComponents(struct node_gra *net, gsl_rng *gen);
void NormalizeCoordinates(struct node_gra *net);
void NormalizeCoordinates3D(struct node_gra *net);
double PotentialEnergy(struct node_gra *net,
		       struct node_box *box_list[], int boxprow);
double PotentialEnergy3D(struct node_gra *net,
			 struct node_box *box_list[], int boxprow);
void CalculateNodeForces(struct node_gra *net, double Fx[], double Fy[],
			 struct node_box *box_list[], int boxprow);
void CalculateNodeForces3D(struct node_gra *net,
			   double Fx[], double Fy[], double Fz[],
			   struct node_box *box_list[], int boxprow);

/*
  -----------------------------------------------------------------------------
  High-level layout functions
  -----------------------------------------------------------------------------
*/
void MDGraphLayout(struct node_gra *net, double drag, double dt,
		   int nsteps, gsl_rng *gen, int nbox);
void MDGraphLayout3D(struct node_gra *net, double drag, double dt,
		     int nsteps, gsl_rng *gen, int nbox);
void MDGraphLayout2Dp(struct node_gra *net, double drag, double dt,
		      int nsteps, gsl_rng *gen, int nbox);

/*
  -----------------------------------------------------------------------------
  Output
  -----------------------------------------------------------------------------
*/
void PrintNodeCoordinates(FILE *outFile, struct node_gra *net);


int CopyCoords(struct node_gra *net_src, struct node_gra *net_des);


#endif /* !RGRAPH_VISUAL_H */
