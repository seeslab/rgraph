/*
  modules.h
  $LastChangedDate$
  $Revision$
*/

#ifndef RGRAPH_MODULES_H
#define RGRAPH_MODULES_H 1

#include <gsl/gsl_rng.h>

/*
  ---------------------------------------------------------------------
  Definition of the group structure
  ---------------------------------------------------------------------
*/
struct group{
  int label;        /* label of the group */
  int size;         /* number of nodes in the group */
  int totlinks;     /* total number of links of the nodes in the group */
  int inlinks;      /* links inside the group */
  int outlinks;     /* links outside the group */
  double totlinksW; /* weighted links of the nodes in the group */
  double inlinksW;  /* weighted links inside the group */
  double outlinksW; /* weighted links outside the group */

  double coorX;
  double coorZ;
  double coorY;

  struct node_lis *nodeList; /* list of nodes in the group */
  struct group *next;        /* next group */

  struct group *offspr;      /* partition of this group */
};


/*
  ---------------------------------------------------------------------
  Group creation and memory allocation
  ---------------------------------------------------------------------
*/
struct group *CreateHeaderGroup();
struct group *CreateGroup(struct group *part, int label);

/*
  ---------------------------------------------------------------------
  Partition creation
  ---------------------------------------------------------------------
*/
struct group *FCreatePartition(FILE *inF);
struct group *CreateEquiNPartition(struct node_gra *net, int gsize);
struct group *CreateEquiNPartitionSoft(int ngroups, int gsize);
struct group *CreatePartitionFromInGroup(struct node_gra *net);

/*
  ---------------------------------------------------------------------
  Partition removal
  ---------------------------------------------------------------------
*/
void RemovePartition(struct group *part);

/*
  ---------------------------------------------------------------------
  Node-group functions
  ---------------------------------------------------------------------
*/
struct node_lis *AddNodeToGroup(struct group *g, struct node_gra *node);
struct node_lis *AddNodeToGroupSoft(struct group *g, char *label);
int RemoveNodeFromGroup(struct group *g, struct node_gra *node);
int MoveNode(struct node_gra *node,
      	     struct group *old,
			 struct group *new);

/* a la DB Stouffer */
struct node_lis *AddNodeToGroupFast(struct group *g, struct node_gra *node);
int RemoveNodeFromGroupFast(struct group *g, struct node_gra *node);
int MoveNodeFast(struct node_gra *node,
		 struct group *old,
		 struct group *new);

/*
  ---------------------------------------------------------------------
  Group and partition operations
  ---------------------------------------------------------------------
*/
struct group *CompressPart(struct group *part);
struct group *GetGroup(struct group *part, int label);
int NGroups(struct group *part);
int NNonEmptyGroups(struct group *part);
int PartitionSize(struct group *part);
void RemoveWithinGroupLinks(struct group *g, int symmetric_sw);
void RemoveBetweenGroupLinks(struct group *part, int symmetric_sw);
void BlockModel(struct group *part,
		    char type_sw,
		    int list_sw);
int NLinksToGroup(struct node_gra* node, struct group *g);
int NWeightLinksToGroup(struct node_gra* node, struct group *g, double w);
int NLinksToGroupByNum(struct node_gra* node, int gLabel);
double StrengthToGroup(struct node_gra* node, struct group *g);
double StrengthToGroupByNum(struct node_gra* node, int gLabel);
int NG2GLinks(struct group *g1, struct group *g2);
int NWeightG2GLinks(struct group *g1, struct group *g2, double w);
double NG2GLinksWeight(struct group *g1, struct group *g2);
void MergeGroups(struct group *g1, struct group *g2);
struct group *CopyGroup(struct group *copy_root, struct group *g);
struct group *CopyPartition(struct group *original);
struct node_gra *BuildNetFromGroup(struct group *group);
struct node_gra *BuildNetFromGroupNeig(struct group *group);
void GroupSizeStatistics(struct group *part,
			 double *theMean,
			 double *theStddev,
			 double *theMin,
			 double *theMax);
struct group *GetEmptyGroup(struct group *part);

/* a la DB Stouffer */
void MergeGroupsFast(struct group *g1, struct group *g2);

/*
  ---------------------------------------------------------------------
  Network-partition operations
  ---------------------------------------------------------------------
*/
void ResetNetGroup(struct node_gra *net);
struct group *ClustersPartition(struct node_gra *net);
void MapPartToNet(struct group *part, struct node_gra *net);
void MapPartToNetSoft(struct group *part, struct node_gra *net);
struct group *ClustersPartition(struct node_gra *net);
void RemoveInterGroupLinks(struct node_gra *net);

/* a la DB Stouffer */
void MapPartToNetFast(struct group *part, struct node_gra *net);

/*
  ---------------------------------------------------------------------
  Group and partition output
  ---------------------------------------------------------------------
*/
void FPrintPartition(FILE *outf, struct group *partition, int list_sw);
void FPrintPajekPartitionFile(char *fname, struct node_gra *net);

/*
  ---------------------------------------------------------------------
  Partition comparison
  ---------------------------------------------------------------------
*/
double MutualInformation(struct group *part1, struct group *part2, int label_sw);
double CorrectlyClassified(struct group *refpart,
			   struct group *actpart);

/*
  ---------------------------------------------------------------------
  Module indentification
  ---------------------------------------------------------------------
*/
double Modularity(struct group *part);
double ModularityWeight(struct group *part);

/*
  ---------------------------------------------------------------------
  Roles
  ---------------------------------------------------------------------
*/
double ParticipationCoefficient(struct node_gra *node);
double WeightedParticipationCoefficient(struct node_gra *node,
										struct group *part);
double WithinModuleRelativeDegree(struct node_gra *node,
								  struct group *group);
double WithinModuleRelativeStrength(struct node_gra *node,
									struct group *group);

struct group *CatalogRoleIdent(struct node_gra *net,
							   struct group *comm);
struct group *CatalogRoleIdentStrength(struct node_gra *net,
			       struct group *comm);


#endif /* !RGRAPH_MODULES_H */
