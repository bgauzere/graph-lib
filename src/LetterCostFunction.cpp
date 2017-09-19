
#include "LetterCostFunction.h"


double LetterDistanceCost::NodeSubstitutionCost(GNode<CMUPoint,double> * n1,GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2)
{
  double dx = n2->attr.x - n1->attr.x;
  double dy = n2->attr.y - n1->attr.y;
  return _alpha * (dx*dx + dy*dy);
}


double LetterDistanceCost::NodeDeletionCost(GNode<CMUPoint,double> * n1,Graph<CMUPoint,double> * g1)
{
  return _alpha * _tnodes;
}


double LetterDistanceCost::NodeInsertionCost(GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g2)
{
  return _alpha * _tnodes;
}



double LetterDistanceCost::EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2)
{
  double d = e1->attr - e2->attr;
  if (d < 0) d *= -1;
  return (1-_alpha) * d; 
}


double LetterDistanceCost::EdgeDeletionCost(GEdge<double> * e1,Graph<CMUPoint,double> * g1)
{
  return (1-_alpha) * _tedges;
}


double LetterDistanceCost::EdgeInsertionCost(GEdge<double> * e2,Graph<CMUPoint,double> * g2)
{
  return (1-_alpha) * _tedges;
}

