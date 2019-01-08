
public class V_E_V implements Comparable<V_E_V> {
Vertex v1;
Vertex v2;
Edge e1;
public V_E_V(Vertex v1, Vertex v2, Edge e1) {
	this.v1=v1;
	this.v2=v2;
	this.e1=e1;
}

public int compareTo(V_E_V vev) {
	// TODO Auto-generated method stub
	 return this.e1._nEdgeCost-vev.e1._nEdgeCost; 
	}
} 

