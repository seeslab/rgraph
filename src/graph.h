/*
  graph.h
*/

#ifndef RGRAPH_GRAPH_H
#define RGRAPH_GRAPH_H 1

#include <search.h>
#include <gsl/gsl_rng.h>

#define MAX_LABEL_LENGTH 100

/*
  ---------------------------------------------------------------------
  Definition of the node_gra structure
  ---------------------------------------------------------------------
*/
struct node_gra{
  char *label;             // label of the node
  int num;                 // ID of the node
  double coorX;
  double coorY;
  double coorZ;
  int state;               // to control wether visited or not
  struct node_lis *neig;   // header of the adjacency list
  struct node_gra *next;   // next node in the graph
  int inGroup;             // group to which the node belongs
  int ivar1;               // number of packets currently at the node
  double dvar1;            // rush variable to store doubles
  int trans;               // translation of the node (for visualization)
  unsigned int degree;     // The length of neig (to speed up modularity computation). 
  double strength;		   // the sum of all link weights (to speed up modularity calculations)
};

/*
  ---------------------------------------------------------------------
  Definition of the node_lis structure
  ---------------------------------------------------------------------
*/
struct node_lis{
  int node;                // num of the referenced node
  char *nodeLabel;         // label of the referenced node
  int status;
  struct node_lis *next;   // next node in the adjacency list
  struct node_gra *ref;    // pointer to the reference node
  double btw;              // link betweenness (rush) variable.
  double weight;           // weight of the link.
};

/*
  ---------------------------------------------------------------------
  Definition of the node_tree structure for the binary tree data
  structure
  ---------------------------------------------------------------------
*/
struct node_tree{
  char *label;                  // label of the node
  struct node_gra *ref;         // pointer to the corresponding node_gra
};

// ---------------------------------------------------------------------
// Structure for the bfs algorithm
// ---------------------------------------------------------------------
struct node_bfs{
  int d;
  struct node_gra *ref;
  struct node_bfs *next;   // Next element in the list
  struct node_bfs *prev;   // Previous element in the list
  struct node_bfs *last;   // Last element in the list
  struct pred *pred;
};

/*
  ---------------------------------------------------------------------
  Predecessor structure for the betweenness calculaton
  ---------------------------------------------------------------------
*/
struct pred{
  struct node_bfs *ref;
  struct pred *next;
};

/*
  ---------------------------------------------------------------------
  Node, link, and graph creation and memory allocation
  ---------------------------------------------------------------------
*/
struct node_gra *CreateHeaderGraph();
struct node_gra *CreateNodeGraph(struct node_gra *net,
				 char *label);
struct node_bfs *CreateHeaderList();
struct node_tree *CreateNodeTree();
int AddAdjacency(struct node_gra *node1,
		 struct node_gra *node2,
		 int auto_link_sw,
		 int add_weight_sw,
		 double weight,
		 int status);
int AddAdjacencySoft(struct node_gra *node1,
		     char *node2Label,
		     int auto_link_sw,
		     int add_weight_sw,
		     double weight,
		     int status);
void RewireAdjacencyByNum(struct node_gra *root);
void RewireAdjacencyByLabel(struct node_gra *root);
void CopyAdjacencyList(struct node_gra *n1,
		       struct node_gra *n2);
struct node_gra *CopyNetwork(struct node_gra *p1);
void *MakeLabelDict(struct node_gra *net);

/*
  ---------------------------------------------------------------------
  Node, link, and graph removal
  ---------------------------------------------------------------------
*/
void FreeNodeTree(struct node_tree *ntree);
void FreeNodeLis(struct node_lis *p);
void FreeAdjacencyList(struct node_lis *p);
void FreeNode(struct node_gra *node);
void RemoveGraph(struct node_gra *p);
void RemoveLink(struct node_gra *n1,
		struct node_gra *n2,
		int symmetric_sw);
void FreeLabelDict(void *dict);

/*
  ---------------------------------------------------------------------
  Network input
  ---------------------------------------------------------------------
*/
int NodeTreeLabelCompare(const void *n1,
			 const void *n2);
struct node_gra *FBuildNetwork(FILE *inFile,
			       int weight_sw,
			       int auto_link_sw,
			       int add_weight_sw,
			       int symmetric_sw);

/*
  ---------------------------------------------------------------------
  Network printing and output
  ---------------------------------------------------------------------
*/
void FPrintDegrees(FILE *file,
		   struct node_gra *p);
void FPrintNetAdjacencyList(FILE *outf,
			    struct node_gra *p,
			    int weight_sw,
			    int symmetric_sw);

void FPrintPajekFile(char *fname,
		     struct node_gra *root,
		     int coor_sw,
		     int weight_sw,
		     int symmetric_sw);

/*
  ---------------------------------------------------------------------
  Node and link operations
  ---------------------------------------------------------------------
*/
struct node_gra *GetNode(int num,
			 struct node_gra *p);
struct node_gra *GetNodeDict(char *label,
			     void *dict);
struct node_lis *GetLink(struct node_gra *n1,
			 int n2);
int IsThereNode(char *label, struct node_gra *p);
int IsThereLink(struct node_gra *n1,
		struct node_gra *n2);
int IsThereLinkSoft(struct node_gra *n1,
		    int n2_num);
int RemoveIsolatedNodes(struct node_gra *root);
void CleanAdjacencies(struct node_gra *net);
void RemoveRandomLinks(struct node_gra *net,
		       int nLinks,
		       int symmetric_sw,
		       gsl_rng *gen);
void AddRandomLinks(struct node_gra *net,
		    int nLinks,
		    int symmetric_sw,
		    gsl_rng *gen);

/*
  ---------------------------------------------------------------------
  BFS list operations
  ---------------------------------------------------------------------
*/
struct node_bfs *GetBFS(struct node_gra *node,
			struct node_bfs *list);
void AddPredecessor(struct node_bfs *node,
		    struct node_bfs *pred);
void ClearPredecessors(struct pred *p);
int CountPredecessors(struct node_bfs *node);
void Enqueue(struct node_gra *node,
	     struct node_bfs *predecessor,
	     struct node_bfs *header,
	     int *size,
	     int dist);
void EnqueueAdjaList(struct node_bfs *lp,
		     struct node_bfs *list,
		     int *size,
		     int d);
struct node_gra *DequeueOne(struct node_bfs *list,
			    struct node_bfs *one,
			    int *size);
struct node_gra *Dequeue(struct node_bfs *list,
			 int *size);
struct node_bfs *RenewQueue(struct node_bfs *list,
			    struct node_bfs *lp,
			    int *size,
			    int d);
void ClearList(struct node_bfs *list,
	       int *size);

/*
  ---------------------------------------------------------------------
  Network resetting
  ---------------------------------------------------------------------
*/
void ResetNodesState(struct node_gra *p);
void RenumberNodes(struct node_gra *net);

/*
  ---------------------------------------------------------------------
  Network randomization
  ---------------------------------------------------------------------
*/
struct node_gra *RandomizeSymmetricNetwork(struct node_gra *net,
					   double times,
					   gsl_rng *gen);

/*
  ---------------------------------------------------------------------
  Network properties
  ---------------------------------------------------------------------
*/
int CountNodes(struct node_gra *p);
unsigned int NodeDegree(struct node_gra *node);
//double SumWeights(struct node_gra *node);
double AverageDegree(struct node_gra *root,
		     int symmetric_sw);
double AverageSquaredDegree(struct node_gra *root);
int TotalNLinks(struct node_gra *p,
		int symmetric_sw);
double NodeStrength(struct node_gra *node);


/* Distance and related functions */
void FPrintDistanceHistogram(FILE *file,
				struct node_gra *root);
void FPrintDistanceHistogramFromNode(FILE *file,
				     struct node_gra *root,
				     int orinode);
double AveragePathLength(struct node_gra *root);
double AverageInverseDistance(struct node_gra *root);

/* Clustering coefficient and related functions */
int SumCommonLinks(struct node_gra *node,
		   struct node_gra *root);
int CalculateLinksBetweenNeig(struct node_bfs *p,
			      struct node_gra *root);
double ClusteringCoefficient(struct node_gra *root);
double ClusteringCoefficient2(struct node_gra *root);
double ClusteringCoefficientDegree(struct node_gra *root);
double SquareClustering(struct node_gra *root);
double OneNodeClusteringCoefficient(struct node_gra *node,
				    struct node_gra *root);
double OneNodeSquareClustering(struct node_gra *node,
			       struct node_gra *root);

/* Betweenness centrality and related functions */
void CalculateLinkBetweenness(struct node_gra *root);
void CalculateBiggestLinkBetweenness(struct node_gra *root,
				     int *n1,
				     int *n2);
void CalculateNodeBetweenness(struct node_gra *net);
void NodeBetweennessStatistics(struct node_gra *net,
			       double *theMean,
			       double *theStddev,
			       double *theMin,
			       double *theMax);

/* Degree correlations */
double Assortativity(struct node_gra *net);
double CalculateKnn(struct node_gra *node);
double CalculateGroupConCor(struct node_gra *net);

/* Components */
int IsGraphConnected(struct node_gra *p);
int AreConnectedList(struct node_gra *root,
		     struct node_gra *n1,
		     int cluslis[]);
int CountStronglyConnectedSets(struct node_gra *root);
struct node_gra *GetLargestStronglyConnectedSet(struct node_gra *root,
						int thres);
struct node_gra *GetLargestWeaklyConnectedSet(struct node_gra *root,
					      int thres);
int GetAllConnectedSets(struct node_gra *network,
                        struct node_gra **net_list); //JORDI
/*
  ---------------------------------------------------------------------
  Node and network comparison
  ---------------------------------------------------------------------
*/
int CommonNeighbors(struct node_gra *n1, struct node_gra *n2);
double JaccardIndex(struct node_gra *n1, struct node_gra *n2);
double TopologicalOverlap(struct node_gra *n1,
			  struct node_gra *n2);
void CompareTwoNetworks(struct node_gra *netA, struct node_gra *netB,
			int *nA_nod, int *nB_nod, int *ncom_nod, 
			double *p_nAB, double *p_nBA,
			int *nA_lin, int *nB_lin, int *ncom_lin, 
			double *p_lAB, double *p_lBA);

/*
  ---------------------------------------------------------------------
  Network operations
  ---------------------------------------------------------------------
*/
struct node_gra * AddTwoNetworks(struct node_gra *netA, struct node_gra *netB);


#endif /* !RGRAPH_GRAPH_H */
