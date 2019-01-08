public class PathSegment {
 protected Vertex _vertex; // the vertex in this path segment
 protected Edge _edge; // the edge associated with this vertex
 public PathSegment(Vertex _vertex,Edge _edge) {
	this._vertex = _vertex; 
	this._edge = _edge; 
 }
 public Vertex getVertex( )
 {
 return _vertex;
 }
 public Edge getEdge( )
 {
 return _edge;
 
 }
 public String toString() {
	 if(_edge != null)
		 return _vertex._strUniqueID.toString() +" " + _edge._nEdgeCost;
	 else
		 return _vertex._strUniqueID.toString() +" NULL";
 }
} 