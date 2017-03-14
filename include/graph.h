#ifndef __PGRAPHH__
#define __PGRAPHH__

#include <vector>
#include <map>

#include <fstream>
#include <cstdlib>
#include <algorithm>    // std::random_shuffle
#include <iostream>

#include <tinyxml.h>
#include "utils.h"

template<class EdgeAttribute>
class GEdge{
private: 
  /** The next neighbourhood node. */
  GEdge *next; 
   /** The rank of the node in the graph array. */
  int incident_node; //GNode ?
  /** The number which identifies the object. */
   int edge_id;

public:
   /** The attribute of the edge */
  EdgeAttribute attr;
  
   /**
    * Creates a new edge to the specified node, with
    * the specified weight and the specified next edge.
    * @param n	the connected node.
    * @param adj the next edge.
    * @param attr the label.
    */
   GEdge( int n, GEdge *adj, EdgeAttribute attr ): next(adj), incident_node(n), edge_id(-1), attr(attr) {};

   /**
    * Creates a new edge to the specified incident_node, with
    * the specified weight and the specified next edge.
    * @param n	the connected incident_node.
    * @param adj	the next edge.
    * @param i       the identifier of the edge
    * @param w	the weight.
    */
   GEdge( int n, GEdge *adj, int i, EdgeAttribute attr): next(adj), incident_node(n), edge_id(i), attr(attr) {};

   /**
    * Deletes the edge
    */
   ~GEdge(){};

   /**
    * Returns the number of the incident incident_node.
    * @return	the number of the connected incident_node.
    */
   int IncidentNode() const { return incident_node; };

  void  setIncidentNode(int new_node)  {  this->incident_node = new_node; };

   /**
    * Returns the next neighbourhood incident_node.
    * @return	the next edge.
    */
  GEdge<EdgeAttribute>* Next() const { return next; };

   /**
    * Sets a new neighbourhood incident_node.
    * @param n	the  neighbourhood edge.
    * @return	the neighbourhood edge.
    */
  GEdge<EdgeAttribute>* Next( GEdge<EdgeAttribute>* n ) { return next=n; };

   /**
    * Returns the index of the referenced object.
    * @return	the index.
    */
   int EdgeId() const { return edge_id; }

   /**
    * Sets the new index of the referenced object.
    * @param i	the new index.
    * @return	the new index of the object.
    */
   int EdgeId( int i ) { return edge_id=i; }
};

/** @brief A node of a graph.
 *
 * The class <code>GNode</code> defines a node.
 * A node indexes an object in a separate array of objects
 * at specified coordinates in the related image.
 * It is characterized by a value and a list
 * of connected nodes.
 */
template <class NodeAttribute, class EdgeAttribute>
class GNode{
private:
  /** The list of the adjacent nodes. */
  GEdge<EdgeAttribute> *adjacents;
  /** The number which identifies the object. */
  int item;

public :
  /** The valuation of the node. */
  NodeAttribute attr;	

  /**
   * Creates a new node with the specified index,
   * and the specified coordinates.
   * @param i	the index of the referenced object.
   * @param p	the specified coordinates.
   */
  GNode<NodeAttribute, EdgeAttribute>( int i, NodeAttribute attr ): adjacents(0), item(i), attr(attr) { };

  /*
   * Node destructor.
   * -> Destroy the list of adjacent node,
   * without worrying about linked nodes.
   * What about attr ?
   */
  ~GNode(){
    GEdge<EdgeAttribute> *q,*p=adjacents;
    
    while ((q=p)) {
      p=p->Next();
      delete q;
    }
  };

  /**
   * Returns the list of all the connected nodes.
   * @return	the list of connected nodes.
   */
  GEdge<EdgeAttribute> * getIncidentEdges() const { return adjacents; };

  GEdge<EdgeAttribute> * Connect( int incidentNode, EdgeAttribute label ) { //XXX: Check for link already here
    return ( adjacents=new GEdge<EdgeAttribute>( incidentNode, adjacents, label ) );
  };
  
  GEdge<EdgeAttribute> * Connect( int incidentNode, int edge_id, EdgeAttribute attr){
    return ( adjacents=new GEdge<EdgeAttribute>( incidentNode, adjacents, edge_id, attr ) );
  };

  int Degree(){
    int degree = 0;
    GEdge<EdgeAttribute>* p = adjacents;
    while(p) {degree ++;p=p->Next();}
    return degree;	       
  };
  /**
   * Deletes the specified node from the list of connected nodes.
   * @param n	the specified node.
   * @return the new list of edges.
   */
  GEdge<EdgeAttribute>* UnConnect( int incidentNode ){
    GEdge<EdgeAttribute> *p = getIncidentEdges();
    GEdge<EdgeAttribute> *q;
   
    if (!p) return NULL;
    if(p->IncidentNode() == incidentNode){
      adjacents = p->Next();
      delete p;
    }
    return adjacents;
   
    while (p->Next())
      {
	if(p->Next()->IncidentNode() == incidentNode){
	  q = p->Next();
	  p->Next(q->Next());
	  delete q;
	}else
	  p = p->Next();
      }
    return adjacents;
  };

  /**
   * Returns the index of the referenced object.
   * @return	the index.
   */
  int Item() const { return item; }

  /**
   * Sets the new index of the referenced object.
   * @param i	the new index.
   * @return	the new index of the object.
   */
  int Item( int i ) { return item=i; }
};

/** @brief A 2D graph.
 *
 * A graph is a set of nodes connected to some other nodes.
 * The two types of graph are supported: directed and undirected;
 * the type must be specified with the constructor.
 * A node is characterized by a value which determines if the
 * node is active in the representation of the graph or not.
 * A node can be any of the Pandore object. It indexes an 
 * item in a separate array of items. (The secret is to use an
 * item as a pointer to a specific objet in an array [no type]).
 * <br>For the use of Graph2d see @ref graph_page.
 */
template< class NodeAttribute, class EdgeAttribute>
class Graph {
private :
  std::vector<GNode<NodeAttribute, EdgeAttribute> *> tnode;
  int nbNodes;
  int nbEdges;
  bool _directed;

  friend class GEdge<EdgeAttribute>;
  
public :
  /**
   * Constructor from a gxl file
   */
  Graph(const char * filename, NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
	EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem));
  
  
  void GraphLoadGXL(const char * filename,
		    NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
		    EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem));
  
  /**
   * Deletes the graph.
   */
  ~Graph(){
    for (int i=0;i<nbNodes;i++) {
      if (tnode[i])
	delete tnode[i];
    }
  };
  
  /**
   * @return true if the graph is directed.
   */
  bool isDirected() const { return _directed; }
  /**
   * Returns the number of nodes.
   * @return	the size.
   */
  int Size() { return nbNodes; };
  int getNbEdges() const { return _directed?nbEdges:nbEdges/2; }
  /**
   * Creates a new graph with no data.
   * @param directed true for creating a directed graph.
   */
  Graph( bool directed =false): tnode(0), nbNodes(0), nbEdges(0), _directed(directed) { }

  
  /**
   * Build a graph according to information encoded in filename.
   * The file corresponding to filename must be in ct (chemDraw) format.
   */
    
  /**
   * Returns the node at the specified coordinates.
   * @param pos	the coordinates.
   * @return	the node at the specified coordinates.
   */
  GNode<NodeAttribute, EdgeAttribute> *operator[]( int pos ){ return(tnode[pos]); }

  /**
   * Returns the node at the specified coordinates.
   * @param pos	the coordinates.
   * @return	the node at the specified coordinates.
   */
  const GNode<NodeAttribute, EdgeAttribute> *operator[]( int pos ) const { return(tnode[pos]); }
   
  /**
   * Adds the specified node that references the specified item. 
   * @param node	the new Node
   * @return	index du noeud.
   */
  int Add( GNode<NodeAttribute, EdgeAttribute>* node ){
    tnode.push_back(node);
    nbNodes ++;
    return tnode.size();
  };
      
  /**
   * Deletes the specified node from the graph. Unlinks it
   * from connected nodes.
   * @param s	the node to be deleted.
   * @return	SUCCESS or FAILURE.
   */
  GNode<NodeAttribute, EdgeAttribute> * Del( int s ){
    GNode<NodeAttribute, EdgeAttribute> *oldNode;
    oldNode = tnode[s];
    tnode[s] = NULL;
    nbNodes --;
    return oldNode;
  };

  GEdge<EdgeAttribute> * Link(int firstNode, int secondNode , EdgeAttribute label){
    GEdge<EdgeAttribute> * e = NULL;
    if (tnode[firstNode] != NULL && tnode[secondNode] != NULL){
      e = tnode[firstNode]->Connect(secondNode,nbEdges,label);
      nbEdges ++;
      if(!_directed){
	tnode[secondNode]->Connect(firstNode,nbEdges,label);
	nbEdges ++;
      }
    }
    return e;
  };
  //TODO : UnLink !

  bool isLinked(int firstNode, int secondNode){
   return (getEdge(firstNode,secondNode) != NULL);      
  };
  
  GEdge<EdgeAttribute> * getEdge(int firstNode, int secondNode){
    GEdge<EdgeAttribute> *p = tnode[firstNode]->getIncidentEdges();
    while(p){
      if (p->IncidentNode() == secondNode)
	return p;
      else
	p=p->Next();
    }
    return NULL;      
  };
  
  GEdge<EdgeAttribute> * getSymmetricEdge(int nodeId,const GEdge<EdgeAttribute> * p){
   if(!_directed){
      return getEdge(p->IncidentNode(), nodeId);
   }
   return NULL;
  };
  void shuffleize(){
    // set some values:
    std::vector<int> perm;
    for (int i=0; i<nbNodes; ++i) perm.push_back(i);
    // using built-in random generator:
    std::random_shuffle ( perm.begin(), perm.end() );
    std::vector<GNode<int,int>*> new_tnode;
    std::vector<int> inv_tnode(Size());
    for(int i=0;i<nbNodes;i++){
      new_tnode.push_back(tnode[perm[i]]);
      inv_tnode[perm[i]] = i;
    }
    for(int i=0;i<nbNodes;i++){
      GEdge<EdgeAttribute> *p = tnode[i]->getIncidentEdges();
      while(p){
	p->setIncidentNode(inv_tnode[p->IncidentNode()]);
	p = p->Next();
      }
    }
    for(int i=0;i<nbNodes;i++)tnode[i] = new_tnode[i];
  }
};

template < class NodeAttribute, class EdgeAttribute>
void Graph<NodeAttribute,EdgeAttribute>::GraphLoadGXL(const char * filename,
						       NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
						       EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem)){
  std::ifstream file(filename,std::ios::in);
  std::vector<char*> v;
  _directed = false;
  TiXmlDocument doc(filename );
  if(!doc.LoadFile()){
    std::cerr << "erreur lors du chargement" << std::endl;
    std::cerr << "error #" << doc.ErrorId() << " : " << doc.ErrorDesc() << std::endl;
  }
   
  TiXmlHandle hdl(&doc);
  std::map<int,int> id_to_index;
  TiXmlElement *elem = hdl.FirstChildElement().FirstChildElement().FirstChildElement().Element();
  while (elem){
    if(strcmp(elem->Value(),"node") == 0){
      int id = std::stoi( elem->Attribute("id"));
      NodeAttribute label = readNodeLabel(elem);
      id_to_index[id] = nbNodes;
      Add(new GNode<NodeAttribute, EdgeAttribute>(id,label));
      
    }else if (strcmp(elem->Value(),"edge") == 0){
        int from=-1;
        int to=-1;

	from = std::stoi(elem->Attribute("from"));
	to = std::stoi(elem->Attribute("to"));
	EdgeAttribute label = readEdgeLabel(elem);
	Link(id_to_index[from], id_to_index[to],label);	      
	
    }
    elem = elem->NextSiblingElement(); // iteration
  }
}

template<class NodeAttribute, class EdgeAttribute>
Graph<NodeAttribute,EdgeAttribute>::Graph(const char * filename, NodeAttribute (*readNodeLabel)(TiXmlElement *elem),EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem)){
  GraphLoadGXL(filename, readNodeLabel,readEdgeLabel);
}


/**
 * Initialize a graph nb_nodes and corresponding to adjacency matrix am
 * the format of am is the following : 
 * am[i][j] : label of edge between nodes i and j, 0 if no edge
 * am[i][i] : label of node i
 */
// static
// 


#endif // __PGRAPHH__
