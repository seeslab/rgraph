// ---------------------------------------------------------------------
// Definition of the group structure
// ---------------------------------------------------------------------
struct group{
  int label;        // label of the group
  int size;         // number of nodes in the group
  int totlinks;     // total number of links of the nodes in the group
  int inlinks;      // links inside the group
  int outlinks;     // links outside the group
  double totlinksW; // wighted links of the nodes in the group
  double inlinksW;  // wighted links inside the group
  double outlinksW; // weighted links outside the group

  double coorX;
  double coorZ;
  double coorY;

  struct node_lis *nodeList; // list of nodes in the group
  struct group *next;       // next group

  struct group *offspr;    // partition of this group
};


// ---------------------------------------------------------------------
// Group creation and memory allocation
// ---------------------------------------------------------------------
struct group *CreateHeaderGroup();
struct group *CreateGroup(struct group *part, int label);

// ---------------------------------------------------------------------
// Partition creation
// ---------------------------------------------------------------------
struct group *FCreatePartition(FILE *inF);
struct group *CreateEquiNPartition(struct node_gra *net, int gsize);
struct group *CreateEquiNPartitionSoft(int ngroups, int gsize);
struct group *CreatePartitionFromInGroup(struct node_gra *net);

// ---------------------------------------------------------------------
// Partition removal
// ---------------------------------------------------------------------
void RemovePartition(struct group *part);

// ---------------------------------------------------------------------
// Node-group functions
// ---------------------------------------------------------------------
struct node_lis *AddNodeToGroup(struct group *g, struct node_gra *node);
struct node_lis *AddNodeToGroupSoft(struct group *g, char *label);
int RemoveNodeFromGroup(struct group *g, struct node_gra *node);
int MoveNode(struct node_gra *node,
	     struct group *old,
	     struct group *new);

// ---------------------------------------------------------------------
// Group and partition operations
// ---------------------------------------------------------------------
struct group *CompressPart(struct group *part);
struct group *GetGroup(struct group *part, int label);
int NGroups(struct group *part);
int NNonEmptyGroups(struct group *part);
int PartitionSize(struct group *part);
void RemoveWithinGroupLinks(struct group *g, int symmetric_sw);
void RemoveBetweenGroupLinks(struct group *part, int symmetric_sw);
double **BlockModel(FILE *outf,
		    struct group *part,
		    char type_sw,
		    int list_sw);
int NLinksToGroup(struct node_gra* node, struct group *g);
double StrengthToGroup(struct node_gra* node, struct group *g);
int NG2GLinks(struct group *g1, struct group *g2);
double NG2GLinksWeight(struct group *g1, struct group *g2);
void MergeGroups(struct group *g1, struct group *g2);
struct group *CopyGroup(struct group *copy_root, struct group *g);
struct group *CopyPartition(struct group *original);
struct node_gra *BuildNetFromPart(struct group *part);

// ---------------------------------------------------------------------
// Network-partition operations
// ---------------------------------------------------------------------
void ResetNetGroup(struct node_gra *net);
struct group *ClustersPartition(struct node_gra *net);
void MapPartToNet(struct group *part, struct node_gra *net);
void MapPartToNetSoft(struct group *part, struct node_gra *net);

// ---------------------------------------------------------------------
// Group and partition output
// ---------------------------------------------------------------------
void FPrintPartition(FILE *outf, struct group *partition, int list_sw);

// ---------------------------------------------------------------------
// Module indentification
// ---------------------------------------------------------------------
double Modularity(struct group *part);
double ModularityWeight(struct group *part);
struct group *SAGroupSplit(struct group *targ,
			   double Ti, double Tf, double Ts,
			   int cluster_sw,
			   struct prng *gen);
