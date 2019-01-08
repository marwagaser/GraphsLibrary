import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.Vector;

public class Graph extends GraphException {

	public ArrayList<ArrayList<Edge>> Matrix = new ArrayList<ArrayList<Edge>>();
	public ArrayList<Vertex> idMatrix = new ArrayList<Vertex>();
	public  ArrayList<String> visited = new ArrayList<String>(); // ,marwa
	public  Vector<PathSegment> dfsPath = new Vector<PathSegment>(); // marwa
	public  Stack<String> tempDfs = new Stack<String>(); // marwa
	public ArrayList<V_E_V> vertexEdgevertex = new ArrayList<V_E_V>(); // contains all the inserted edges and their
																		// vertices
	public ArrayList<V_E_V> MSTvev = new ArrayList<V_E_V>(); // temp result
	public ArrayList<Integer> costList = new ArrayList<Integer>();// cost of each edge in the vertexEdgevertex
	public Vector<PathSegment> minSpan = new Vector<PathSegment>(); // vector to store min spanning tree

	class Answer {
		public Double min;
		Vertex[] v;

		Answer() {
			min = Double.MAX_VALUE;
			v = new Vertex[2];
		}
	}

	public double pythagoreanDistance(double x1, double y1, double x2, double y2) {
		return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
	}

	public Graph(String strMessage) {
		super(strMessage);
		// TODO Auto-generated constructor stub
	}

	public int getIndex(String id) throws GraphException {
		for (int i = 0; i < idMatrix.size(); i++) {
			if ((idMatrix.get(i)._strUniqueID.toString()).equals(id)) {
				return i;
			}
		}
		throw new GraphException("Doesn't exist");
	}

	public String getLibraryName() {
		return "UndirectedGraphs";
	}

	public String getLibraryVersion() {
		return "0.1";
	}

	public void insertVertex(String strUniqueID, String strData, int nX, int nY) throws GraphException {
		if(strUniqueID == null)
			throw new GraphException("Unexpected null in input");
		StringBuffer strId = new StringBuffer();
		strId.append(strUniqueID);
		StringBuffer strDt = new StringBuffer();
		strDt.append(strData);
		for (int i = 0; i < idMatrix.size(); i++) {
			Matrix.get(i).add(new Edge());

		}
		Vertex x = new Vertex(strId, strDt, nX, nY);
		idMatrix.add(x);
		visited.add("unexplored"); // marwa
		ArrayList<Edge> list = new ArrayList<Edge>();

		for (int i = 0; i < idMatrix.size(); i++) {
			list.add(new Edge());
		}

		Matrix.add(list);

	}

	public void insertEdge(String strVertex1UniqueID, String strVertex2UniqueID, String strEdgeUniqueID,
			String strEdgeData, int nEdgeCost) throws GraphException {
		if(strVertex1UniqueID == null || strVertex2UniqueID == null || strEdgeUniqueID == null)
			throw new GraphException("Unexpected null in input");
		boolean flag = false;
		int firstvertex = -2;
		int secondvertex = -1;
		firstvertex = getIndex(strVertex1UniqueID);
		secondvertex = getIndex(strVertex2UniqueID);
		StringBuffer ID = new StringBuffer();
		ID.append(strEdgeUniqueID);

		StringBuffer Data = new StringBuffer();
		Data.append(strEdgeData);

		Edge e = new Edge(ID, Data, nEdgeCost);
		for (int i = 0; i < 2; i++) {
			Matrix.get(firstvertex).get(secondvertex)._nEdgeCost = nEdgeCost;
			Matrix.get(firstvertex).get(secondvertex)._strData = Data;
			Matrix.get(firstvertex).get(secondvertex)._strUniqueID = ID;
			int temp = firstvertex;
			firstvertex = secondvertex;
			secondvertex = temp;
		}
		if (flag == false) {
			Vertex v1 = findVertex(strVertex1UniqueID);// marwa
			Vertex v2 = findVertex(strVertex2UniqueID);// marwa
			Edge e1 = findEdge(strEdgeUniqueID); // marwa
			V_E_V vev = new V_E_V(v1, v2, e1);// marwa
			vertexEdgevertex.add(vev);
			costList.add(e1._nEdgeCost);
			flag = true;
		}
	}

	//////////////////// MARWA ///////////////////////
	public Vertex findVertex(String vertexId) throws GraphException { // method to find vertex given id
		if(vertexId == null)
			throw new GraphException("Unexpected null in input");
		boolean flag = false;
		for (int i = 0; i < idMatrix.size(); i++) {
			if (idMatrix.get(i)._strUniqueID.toString().equals(vertexId)) {
				flag = true;
				return idMatrix.get(i);
			}
		}
		if (flag == false) {
			throw new GraphException("no such vertex");
		}
		return null;
	}

	public Edge findEdge(String strEdgeUniqueID) throws GraphException { // method to find edge given id
		if(strEdgeUniqueID == null)
			throw new GraphException("Unexpected null in input");
		boolean flag = false;
		for (int i = 0; i < Matrix.size(); i++) {

			for (int j = 0; j < Matrix.get(0).size(); j++) {
				if ((Matrix.get(i).get(j)._strUniqueID != null)
						&& Matrix.get(i).get(j)._strUniqueID.toString().equals(strEdgeUniqueID)) {
					flag = true;
					// System.out.println( Matrix.get(i).get(j)._strUniqueID.toString());
					return Matrix.get(i).get(j);

				}
			}

		}
		if (!flag)
			throw new GraphException("This edge does not exist!");
		return null;
	}
	//////////////// MARWA

	public void removeVertex(String strVertexUniqueID) throws GraphException {
		if(strVertexUniqueID == null)
			throw new GraphException("Unexpected null in input");
		int i;
		boolean flag = false;
		for (i = 0; i < idMatrix.size(); i++) {
			if (idMatrix.get(i)._strUniqueID.toString().equals(strVertexUniqueID)) {
				idMatrix.remove(i);
				flag = true;
				break;
			}
		}
		i--;
		Matrix.remove(i);

		for (ArrayList<Edge> row : Matrix) {
			row.remove(i);
		}
		if (!flag)
			throw new GraphException("This vertex does not exist!");
	}

	public void removeEdge(String strEdgeUniqueID) throws GraphException {
		if(strEdgeUniqueID == null)
			throw new GraphException("Unexpected null in input");
		boolean flag = false;
		for (int i = 0; i < Matrix.size(); i++) {

			for (int j = 0; j < Matrix.get(0).size(); j++) {
				if ((Matrix.get(i).get(j)._strUniqueID != null)
						&& Matrix.get(i).get(j)._strUniqueID.toString().equals(strEdgeUniqueID)) {
					Matrix.get(i).get(j)._nEdgeCost = 0;
					Matrix.get(i).get(j)._strData = null;
					Matrix.get(i).get(j)._strUniqueID = null;
					flag = true;
				}
			}

		}
		if (!flag)
			throw new GraphException("This edge does not exist!");

	}

	public Vector<Edge> incidentEdges(String strVertexUniqueID) throws GraphException {
		if(strVertexUniqueID == null)
			throw new GraphException("Unexpected null in input");
		Vector<Edge> vector = new Vector<Edge>();
		int index = getIndex(strVertexUniqueID);
		for (int i = 0; i < idMatrix.size(); i++) {
			Edge e1 = Matrix.get(i).get(index);
			Edge e2 = Matrix.get(index).get(i);
			if (e1._strUniqueID != null)
				vector.add(Matrix.get(i).get(index));
			if (e2._strUniqueID != null)
				vector.add(Matrix.get(index).get(i));
		}
		return vector;
	}

	public Vector<Vertex> vertices() throws GraphException {
		Vector<Vertex> vector = new Vector<Vertex>();
		for (int i = 0; i < idMatrix.size(); i++) {
			vector.add(idMatrix.get(i));
		}
		return vector;
	}

	public Vector<Edge> edges() throws GraphException {
		Vector<Edge> vector = new Vector<Edge>();
		for (int i = 0; i < idMatrix.size(); i++) {
			for (int j = 0; j < idMatrix.size(); j++) {
				if (Matrix.get(i).get(j)._strUniqueID != null) {
					vector.add(Matrix.get(i).get(j));
				}
			}
		}
		return vector;
	}

	public Vertex[] endVertices(String strEdgeUniqueID) throws GraphException {
		if(strEdgeUniqueID == null)
			throw new GraphException("Unexpected null in input");
		for (int i = 0; i < idMatrix.size(); i++) {
			for (int j = 0; j < idMatrix.size(); j++) {
				if (Matrix.get(i).get(j)._strUniqueID != null
						&& Matrix.get(i).get(j)._strUniqueID.toString().equals(strEdgeUniqueID)) {
					Vertex[] v = new Vertex[2];
					v[0] = idMatrix.get(i);
					v[1] = idMatrix.get(j);
					return v;
				}
			}
		}
		throw new GraphException("Does not exist!");
	}

	public Vertex opposite(String strVertexUniqueID, String strEdgeUniqueID) throws GraphException {
		if(strVertexUniqueID == null || strEdgeUniqueID == null)
			throw new GraphException("Unexpected null in input");
		Vertex[] v = endVertices(strEdgeUniqueID);
		if (!(v[0]._strUniqueID.toString().equals(strVertexUniqueID)))
			return v[0];
		else
			return v[1];
	}

	public void printVertices() {
		for (Vertex vertex : idMatrix) {
			System.out.print(vertex._strUniqueID.toString() + " ");
		}
		System.out.println("");
	}

	public void printMatrix() {
		for (ArrayList<Edge> row : Matrix) {
			for (Edge edge : row) {
				if (edge._strUniqueID != null)
					System.out.print(edge._strUniqueID.toString() + " ");
				else
					System.out.print(null + " ");
			}
			System.out.println("");
		}
	}

	public Vertex[] closestPair() throws GraphException {
		// ArrayList<ArrayList<Edge>> matrix1 =new ArrayList<ArrayList<Edge>>(
		// Matrix.subList(0, Matrix.size()/2));
		Answer a = new Answer();
		return closestPairHelper(new ArrayList<Vertex>(idMatrix.subList(0, idMatrix.size() / 2)),
				new ArrayList<Vertex>(idMatrix.subList(idMatrix.size() / 2, idMatrix.size())), a).v;

	}

	public Answer closestPairHelper(ArrayList<Vertex> list1, ArrayList<Vertex> list2, Answer p) throws GraphException {

		p = closestPairCrossing(list1, list2);
		if (list1.size() != 0) {
			Answer a = closestPairHelper(new ArrayList<Vertex>(list1.subList(0, list1.size() / 2)),
					new ArrayList<Vertex>(list1.subList(list1.size() / 2, list1.size())), p);
			if (list2.size() != 0) {
				Answer a2 = closestPairHelper(new ArrayList<Vertex>(list2.subList(0, list2.size() / 2)),
						new ArrayList<Vertex>(list2.subList(list2.size() / 2, list2.size())), p);
				if (a2.min < a.min)
					a = a2;
			}
			if (a.min < p.min) {
				return a;
			} else {
				return p;
			}
		}
		return p;
	}

	public Answer closestPairCrossing(ArrayList<Vertex> l1, ArrayList<Vertex> l2) throws GraphException {
		Answer a = new Answer();
		for (int i = 0; i < l1.size(); i++) {
			for (int j = 0; j < l2.size(); j++) {
				Double num = pythagoreanDistance(l1.get(i)._nX, l1.get(i)._nY, l2.get(j)._nX, l2.get(j)._nY);
				if (num < a.min) {
					a.min = num;
					a.v[0] = l1.get(i);
					a.v[1] = l2.get(j);
				}
			}
		}
		return a;
	}

	public void dfs(String strStartVertexUniqueID, Visitor visitor) throws GraphException {
		if(strStartVertexUniqueID == null || visitor == null)
			throw new GraphException("Unexpected null in input");
		Stack<StringBuffer> s = new Stack<StringBuffer>();
		ArrayList<StringBuffer> tr = new ArrayList<StringBuffer>();
		Vertex vs = this.idMatrix.get(this.getIndex(strStartVertexUniqueID));
		visitor.visit(vs);
		boolean flag = false;
		s.push(vs.getUniqueID());
		while (s.size() != this.idMatrix.size()) {
			// System.out.print(s.peek());
			Vector<Edge> e = this.incidentEdges(vs.getUniqueID().toString());
			Vertex v = this.opposite(vs.getUniqueID().toString(), e.get(0).getUniqueID().toString());
			int i = 0;
			while (s.contains(v.getUniqueID()) && i < e.size()) {
				v = this.opposite(vs.getUniqueID().toString(), e.get(i).getUniqueID().toString());
				// System.out.println(v.getUniqueID().toString());

				i++;
			}
			if (s.contains(v.getUniqueID())) {
				flag = true;
			}
			if (flag == false) {
				s.push(v.getUniqueID());
				visitor.visit(e.get(i));
				visitor.visit(v);

			}
			vs = v;
			if (flag == true) {
				break;
			}
		}
		int o = 1;
		int j = 1;
		boolean f = false;
		Vector<Edge> e;
		Vertex v;
		while (s.size() != this.idMatrix.size()) {
			while (j > 0) {
				tr.add(s.pop());
				// System.out.println(s.peek());
				e = this.incidentEdges(s.peek().toString());
				v = this.opposite(s.peek().toString(), e.get(0).getUniqueID().toString());
				int i = 0;
				while ((s.contains(v.getUniqueID()) || tr.contains(v.getUniqueID())) && (i < e.size())) {
					v = this.opposite(s.peek().toString(), e.get(i).getUniqueID().toString());
					// System.out.println(v.getUniqueID().toString());

					i++;
				}
				if (s.contains(v.getUniqueID()) || tr.contains(v.getUniqueID())) {
					f = true;
					// visitor.visit(e.get(i));
					// visitor.visit(v);
				}
				if (f == false) {
					s.push(v.getUniqueID());
					// System.out.print(true);
					visitor.visit(e.get(i));
					visitor.visit(v);

				}
				j--;
				o++;
				f = false;
			}
			//
			j = o;
			for (int k = tr.size() - 1; k >= 0; k--) {
				s.push(tr.get(k));
			}
			tr.clear();
		}

	}

	public void bfs(String strStartVertexUniqueID, Visitor visitor) throws GraphException, InterruptedException {
		if(strStartVertexUniqueID == null || visitor == null)
			throw new GraphException("Unexpected null in input");
		ArrayList<String> visitedvertix = new ArrayList<String>();
		ArrayList<String> visitededge = new ArrayList<String>();

		int indexofVertix = getIndex(strStartVertexUniqueID);
		Vertex vertex = idMatrix.get(indexofVertix);
		visitor.visit(vertex);
		visitedvertix.add(vertex._strUniqueID.toString());

		ArrayList<String> bfssorted = new ArrayList<String>();
		LinkedList<Vertex> vertexlist = new LinkedList<Vertex>();

		vertexlist.addLast(vertex);
		bfssorted.add(vertex._strUniqueID.toString());

		while (vertexlist.size() != 0) {
			vertex = vertexlist.removeFirst();

			Vector<Edge> incidentEdges = incidentEdges(vertex._strUniqueID.toString());
			while (!incidentEdges.isEmpty()) {
				Edge edge = incidentEdges.remove(0);
				Vertex vrtx = opposite(vertex._strUniqueID.toString(), edge._strUniqueID.toString());
				boolean f1 = false;
				for (int io = 0; io < visitedvertix.size(); io++) {
					if (visitedvertix.get(io).equals(vrtx._strUniqueID.toString())) {
						f1 = true;
					}
				}
				if (!f1) {
					visitedvertix.add(vrtx._strUniqueID.toString());
					visitededge.add(edge._strUniqueID.toString());
					vertexlist.add(vrtx);
					visitor.visit(edge);
					visitor.visit(vrtx);

				} else {
					f1 = false;
					for (int io = 0; io < visitededge.size(); io++) {
						if (visitededge.get(io).equals(edge._strUniqueID.toString())) {
							f1 = true;
						}
					}
					if (!f1) {
						visitededge.add(edge._strUniqueID.toString());
						visitor.visit(edge);
					}
				}
			}
		}

		for (int k = 0; k < visitedvertix.size(); k++) {
			System.out.println(visitedvertix.get(k));
		}
	}

	public Vector<PathSegment> pathDFS(String strStartVertexUniqueID, String strEndVertexUniqueID)
			throws GraphException {
		if(strStartVertexUniqueID == null || strEndVertexUniqueID == null)
			throw new GraphException("Unexpected null in input");
		int startIndex = getIndex(strStartVertexUniqueID); // get its index
		visited.set(startIndex, "visited");
		tempDfs.push(strStartVertexUniqueID); // push your start vertex
		if (strStartVertexUniqueID.equals(strEndVertexUniqueID)) {
			for (int i = 0; i < tempDfs.size(); i = i + 2) {
				if (i == tempDfs.size() - 1) {
					// PathSegment ps = new PathSegment (findVertex(tempDfs.elementAt(i)),null);
					//dfsPath.add(new PathSegment(findVertex(tempDfs.elementAt(i)), null));
                    //System.out.println("VERTEX: " + findVertex(tempDfs.elementAt(i))._strUniqueID + " EDGE: " + null);
				} else {
					// PathSegment ps = new PathSegment
					// (findVertex(tempDfs.elementAt(i)),findEdge(tempDfs.elementAt(i+1)));
					dfsPath.add(new PathSegment(findVertex(tempDfs.elementAt(i)), findEdge(tempDfs.elementAt(i + 1))));
					System.out.println("VERTEX: " + findVertex(tempDfs.elementAt(i))._strUniqueID + " EDGE: "
							+ findEdge(tempDfs.elementAt(i + 1))._strUniqueID);
				}
			}
			
			return dfsPath;
		}

		for (int i = 0; i < incidentEdges(strStartVertexUniqueID).size(); i++) {
			if (incidentEdges(strStartVertexUniqueID).get(i)._EdgeStatus.equals("unexplored")) {
				// opposite(strStartVertexUniqueID,incidentEdges(strStartVertexUniqueID).get(i)._strUniqueID.toString());
				int oppIndex = getIndex(opposite(strStartVertexUniqueID,
						incidentEdges(strStartVertexUniqueID).get(i)._strUniqueID.toString())._strUniqueID.toString());
				if (visited.get(oppIndex).equals("unexplored")) {
					incidentEdges(strStartVertexUniqueID).get(i)._EdgeStatus = "discovery";
					tempDfs.push(incidentEdges(strStartVertexUniqueID).get(i)._strUniqueID.toString());
					pathDFS(opposite(strStartVertexUniqueID,
							incidentEdges(strStartVertexUniqueID).get(i)._strUniqueID.toString())._strUniqueID
									.toString(),
							strEndVertexUniqueID);
					tempDfs.pop();
				} else {
					incidentEdges(strStartVertexUniqueID).get(i)._EdgeStatus = "back";
				}
			}
		}
		tempDfs.pop();

		return null;
	}

	public Vector<Vector<PathSegment>> findShortestPathBF(String strStartVertexUniqueID) throws GraphException {
		if(strStartVertexUniqueID == null)
			throw new GraphException("Unexpected null in input");
		int[] dist = new int[idMatrix.size()];
		Vector<Vector<PathSegment>> r = new Vector<Vector<PathSegment>>();
		for (int i = 0; i < dist.length; i++) {
			dist[i] = Integer.MAX_VALUE;
			r.add(null);
		}

		int vertexIndex = getIndex(strStartVertexUniqueID);
		dist[vertexIndex] = 0;
		r.remove(vertexIndex);
		r.add(vertexIndex, new Vector<PathSegment>());
		Vector<Edge> edges = edges();
		Vector<Vertex> vertices = vertices();
//		for(int i=0; i< edges.size();i++) {
//		if(endVertices(edges.get(i)._strUniqueID.toString())[0]._strUniqueID.toString().equals(strStartVertexUniqueID) || endVertices(edges.get(i)._strUniqueID.toString())[1]._strUniqueID.toString().equals(strStartVertexUniqueID)) {
//			edges.add(0, edges.remove(i));
//		}
//	}

		for (int i = 1; i < vertices.size(); ++i) {
			for (int j = 0; j < edges.size(); ++j) {
				Vertex u = endVertices(edges.get(j)._strUniqueID.toString())[0];
				Vertex v = endVertices(edges.get(j)._strUniqueID.toString())[1];
				int uIndex = vertices.size();
				int vIndex = vertices.size();
				for (int a = 0; a < vertices.size(); a++) {
					if (vertices.get(a).equals(u))
						uIndex = a;
					if (vertices.get(a).equals(v))
						vIndex = a;
				}
				if (uIndex == vertices.size() || vIndex == vertices.size()) {
					throw new GraphException("Could not find vertex index");
				}
				int weight = edges.get(j).getCost();
				
				if (uIndex == vertexIndex && weight < dist[vIndex]) {
					dist[vIndex] = weight;
					r.remove(vIndex);
					r.add(vIndex, new Vector<PathSegment>());
					r.get(vIndex).add(new PathSegment(u, edges.get(j)));
				} else if (vIndex == vertexIndex && weight < dist[uIndex]) {
					dist[uIndex] = weight;
					r.remove(uIndex);
					r.add(uIndex, new Vector<PathSegment>());
					r.get(uIndex).add(new PathSegment(v, edges.get(j)));
				} else if (dist[uIndex] != Integer.MAX_VALUE && weight + dist[uIndex] < dist[vIndex]) {
					dist[vIndex] = weight + dist[uIndex];
					Vector<PathSegment> n = new Vector<PathSegment>();
					for (PathSegment pathSegment : r.get(uIndex)) {
						n.add(pathSegment);
					}
					n.add(new PathSegment(u, edges.get(j)));
					r.remove(vIndex);
					r.add(vIndex, n);
				} else if (dist[vIndex] != Integer.MAX_VALUE && weight + dist[vIndex] < dist[vIndex]) {
					dist[uIndex] = weight + dist[vIndex];
					Vector<PathSegment> n = new Vector<PathSegment>();
					for (PathSegment pathSegment : r.get(vIndex)) {
						n.add(pathSegment);
					}
					n.add(new PathSegment(v, edges.get(j)));
					r.remove(uIndex);
					r.add(uIndex, n);
				}
				
				int tempI = uIndex;
				uIndex = vIndex;
				vIndex = tempI;
				
				Vertex tempV = u;
				u = v;
				v = tempV;
				if (uIndex == vertexIndex && weight < dist[vIndex]) {
					dist[vIndex] = weight;
					r.remove(vIndex);
					r.add(vIndex, new Vector<PathSegment>());
					r.get(vIndex).add(new PathSegment(u, edges.get(j)));
				} else if (vIndex == vertexIndex && weight < dist[uIndex]) {
					dist[uIndex] = weight;
					r.remove(uIndex);
					r.add(uIndex, new Vector<PathSegment>());
					r.get(uIndex).add(new PathSegment(v, edges.get(j)));
				} else if (dist[uIndex] != Integer.MAX_VALUE && weight + dist[uIndex] < dist[vIndex]) {
					dist[vIndex] = weight + dist[uIndex];
					Vector<PathSegment> n = new Vector<PathSegment>();
					for (PathSegment pathSegment : r.get(uIndex)) {
						n.add(pathSegment);
					}
					n.add(new PathSegment(u, edges.get(j)));
					r.remove(vIndex);
					r.add(vIndex, n);
				} else if (dist[vIndex] != Integer.MAX_VALUE && weight + dist[vIndex] < dist[vIndex]) {
					dist[uIndex] = weight + dist[vIndex];
					Vector<PathSegment> n = new Vector<PathSegment>();
					for (PathSegment pathSegment : r.get(vIndex)) {
						n.add(pathSegment);
					}
					n.add(new PathSegment(v, edges.get(j)));
					r.remove(uIndex);
					r.add(uIndex, n);
				}
			}
		}
		for (int i = 0; i < dist.length; i++) {
			System.out.println(vertices.get(i)._strUniqueID.toString() + " " + dist[i]);
		}
		return r;
	}

	public Vector<Vector<PathSegment>> findAllShortestPathsFW() throws GraphException {
		int dist[][] = new int[idMatrix.size()][idMatrix.size()];

		Vector<Vector<Vector<PathSegment>>> fw3d = new Vector<Vector<Vector<PathSegment>>>();

		/*
		 * Initialize the solution matrix same as input graph matrix. Or we can say the
		 * initial values of shortest distances are based on shortest paths considering
		 * no intermediate vertex.
		 */

		Vector<Vertex> vertices = vertices();
		for (int i = 0; i < vertices.size(); i++) {
			fw3d.add(new Vector<Vector<PathSegment>>());
			for (int j = 0; j < vertices.size(); j++) {
				fw3d.get(i).add(new Vector<PathSegment>());
				Edge e = Matrix.get(i).get(j);
				dist[i][j] = e._strUniqueID == null ? Integer.MAX_VALUE / 2 : e._nEdgeCost;
				/*
				 * for (int k = 0; k< vertices.size(); k++) { Vector<PathSegment> m =
				 * r.get(i).get(j); m.add(null); }
				 */
			}
		}

		/*
		 * Add all vertices one by one to the set of intermediate vertices. ---> Before
		 * start of an iteration, we have shortest distances between all pairs of
		 * vertices such that the shortest distances consider only the vertices in set
		 * {0, 1, 2, .. k-1} as intermediate vertices. ----> After the end of an
		 * iteration, vertex no. k is added to the set of intermediate vertices and the
		 * set becomes {0, 1, 2, .. k}
		 */
		for (int k = 0; k < vertices.size(); k++) {
			// Pick all vertices as source one by one
			for (int i = 0; i < vertices.size(); i++) {
				// Pick all vertices as destination for the
				// above picked source
				if (k == i) {
					dist[k][i] = 0;
					fw3d.get(k).remove(i);
					fw3d.get(k).add(i, new Vector<PathSegment>());
				}
				for (int j = 0; j < vertices.size(); j++) {
					// If vertex k is on the shortest path from
					// i to j, then update the value of dist[i][j]

					if (dist[i][k] + dist[k][j] < dist[i][j]) {
						dist[i][j] = dist[i][k] + dist[k][j];
						Vector<PathSegment> iToK = fw3d.get(i).get(k);
						Vector<PathSegment> kToJ = fw3d.get(k).get(j);
						if (iToK.size() == 0)
							iToK.add(new PathSegment(vertices.get(i), Matrix.get(i).get(k)));
						if (kToJ.size() == 0)
							kToJ.add(new PathSegment(vertices.get(k), Matrix.get(k).get(j)));
						// System.out.println(iToK);
						fw3d.get(i).remove(j);
						fw3d.get(i).add(j, new Vector<PathSegment>());
						for (PathSegment pathSegment : iToK) {
							fw3d.get(i).get(j).add(pathSegment);
						}
						// int index = r.get(i).get(j).size();
						for (PathSegment pathSegment : kToJ) {
							fw3d.get(i).get(j).add(pathSegment);
						}
						// PathSegment link = new PathSegment(vertices.get(k), Matrix[k][]);

						// System.out.println(dist[i][j]);
					}

				}
			}
		}

		// Print the shortest distance matrix
		Vector<Vector<PathSegment>> r = new Vector<Vector<PathSegment>>();
		for (int i = 0; i < fw3d.size(); i++) {
			for (int j = 0; j < fw3d.get(i).size(); j++) {
				// fw3d.get(i).get(j).add(0, new PathSegment(vertices.get(i), null));
				r.add(fw3d.get(i).get(j));

			}
		}
		printSolution(dist);
		return r;
	}

	static void printSolution(int dist[][]) {

		for (int i = 0; i < dist.length; i++) {
			for (int j = 0; j < dist[i].length; j++) {
				System.out.print(dist[i][j] + "   ");
			}
			System.out.println();
		}
	}

	public Vector<PathSegment> minSpanningTree() throws GraphException {
		int vMinus1 = 0;
		TesterGraph tg = new TesterGraph(idMatrix.size());
		while (vMinus1 < idMatrix.size() - 1) {

			int minIndex = costList.indexOf(Collections.min(costList)); // gets minimum cost
			V_E_V vev = vertexEdgevertex.get(minIndex);
			String edge = vev.e1._strUniqueID.toString();
			String vertex1 = vev.v1._strUniqueID.toString();
			String vertex2 = vev.v2._strUniqueID.toString();

			Edge edgeE = vev.e1;
			Vertex vertexV1 = vev.v1;
			Vertex vertexV2 = vev.v2;
			boolean c = tg.addEdge(vertex1, vertex2, tg);
			// System.out.println(c);
			if (c == true) {
				costList.remove(minIndex);
				vertexEdgevertex.remove(minIndex);
			} else {
				vMinus1 += 1;
				MSTvev.add(vev);
				costList.remove(minIndex);
				vertexEdgevertex.remove(minIndex);
			}
		}
		System.out.println("Our res");
		for (int i = 0; i < MSTvev.size(); i++) {
			// System.out.println(MSTvev.get(i).v1._strUniqueID+"
			// "+MSTvev.get(i).v2._strUniqueID +" "+ MSTvev.get(i).e1._nEdgeCost);
			minSpan.add(new PathSegment(null, MSTvev.get(i).e1));
			// System.out.println(minSpan.get(index));
		}
		for (int i = 0; i < minSpan.size(); i++)
			System.out.print(minSpan.get(i)._edge._strUniqueID + " ");
		// tg.isCyclic();
		return minSpan;
	}

	public static void main(String[] args) throws GraphException, InterruptedException {
		/*
		 * Graph g= new Graph("First Graph"); g.insertVertex("a", "1", 0, 0);
		 * g.insertVertex("b", "2", 5, 5); g.insertVertex("e", "2", 2, 5);
		 * g.insertVertex("c","3",20,20); g.insertEdge("a", "b", "first", "first", 1);
		 * g.insertEdge("b", "c", "second", "second", 1); //g.insertEdge("c", "c",
		 * "third", "third", 1); g.insertEdge("c", "a", "fourth", "fourth", 1);
		 * //g.insertEdge("c", "a", "fourth", "fourth", 1); g.insertEdge("a", "e",
		 * "fifth", "fifth", 1); GradingVisitor t= new GradingVisitor(); g.dfs("a",t);
		 * //System.out.print(g.incidentEdges("b").get(2).getData());
		 * //g.printVertices(); System.out.println(t.getResult());
		 */
		// g.printMatrix();
		// System.out.println(g.closestPair()[0]._strUniqueID.toString() + " " +
		// g.closestPair()[1]._strUniqueID.toString());
		// bfs testing
		/*
		 * Graph g= new Graph("First Graph"); g.insertVertex("a", "0", 0, 0);
		 * g.insertVertex("b", "1", 0, 1); g.insertVertex("c", "2", 0, 0);
		 * g.insertVertex("d", "3", 0, 1); g.insertVertex("e", "4", 0, 1);
		 * 
		 * g.insertEdge("c", "d", "3", "3", 1); g.insertEdge("a", "b", "0", "0", 1);
		 * g.insertEdge("b", "d", "2", "2", 1); g.insertEdge("a", "c", "1", "1", 1);
		 * g.insertEdge("d", "e", "4", "4", 1);
		 * 
		 * 
		 * GradingVisitor visitor=new GradingVisitor(); g.bfs("a",visitor);
		 * System.out.println(visitor.getResult());
		 */
		//
		// Graph g = new Graph("First Graph");
		// g.insertVertex("a", "1", 0, 0);
		// g.insertVertex("b", "2", 0, 1);
		// g.insertVertex("c", "3", 1, 0);
		// g.insertVertex("d", "4", 1, 0);
		// g.insertEdge("a", "b", "first", "first", 5);
		// g.insertEdge("b", "c", "second", "second", 3);
		// g.insertEdge("c", "d", "sixth", "sixth", 1);
		// g.insertEdge("a", "d", "fourth", "fourth", 10);
		// System.out.println(g.findAllShortestPathsFW());
		//
		//
		//

		Graph g8= new Graph("First Graph");
//		g8.insertVertex("0", "0", 400, 1);
//		g8.insertVertex("1", "1", 1501, 2);
//		g8.insertVertex("2", "2", 200, 3);
//		g8.insertVertex("3", "3", 1500, 4);
//		g8.insertVertex("4", "4", 500, 0);
//		g8.insertVertex("5", "5", 8, 601);
//		g8.insertVertex("6", "6", 90, 4);
//		g8.insertVertex("7", "7", 2, 100);
//		
//		g8.insertEdge("0", "1", "2", "2", 2);
//		g8.insertEdge("1", "2", "3", "3", 3);
//		g8.insertEdge("2", "3", "4", "4", 4);
//		g8.insertEdge("4", "5", "30", "3", 3);
//		g8.insertEdge("5", "6", "1", "1", 1);
//		g8.insertEdge("6", "7", "31", "3", 3);
//		g8.insertEdge("0", "4", "32", "3", 3);
//		g8.insertEdge("1", "5", "20", "2", 2);
//		g8.insertEdge("2", "6", "21", "2", 2);
//		g8.insertEdge("3", "7", "22", "2", 2);
//		g8.insertEdge("1", "4", "23", "2", 2);
//		g8.insertEdge("3", "6", "24", "2", 2);
		g8.insertVertex("a", "a", 0, 0);
		g8.insertVertex("b", "b", 0, 0);
		g8.insertVertex("c", "c", 0, 0);
		g8.insertVertex("d", "d", 0, 0);
		g8.insertEdge("a", "b", "1", "strEdgeData", 5);
		g8.insertEdge("b", "c", "2", "strEdgeData", 3);
		g8.insertEdge("c", "d", "3", "strEdgeData", 1);
		g8.insertEdge("a", "d", "4", "strEdgeData", 10);
		
		System.out.println(g8.findAllShortestPathsFW());
		// pathDFS("a","b");

		/*
		 * //Test 1 testing opposite working in both directions //The next 4 lines
		 * uncomment only one vertex x and comment the rest
		 * 
		 * //Vertex x = opposite("b","first"); //Vertex x = opposite("b","first");
		 * //Vertex x = opposite("b","second"); Vertex x = opposite("c","second");
		 * System.out.println(x._strUniqueID);
		 */

		/*
		 * //Test 2 Override test g.insertEdge("a", "b", "fourth", "fourth", 1);
		 * //Vertex x = opposite("a","fourth"); //should work since it overrides first
		 * //Vertex x = opposite("a","first"); //should throw an error because it is
		 * overridden by fourth System.out.println(x._strUniqueID);
		 */

		/*
		 * //Test 3 Deleting Vertex test (Uncomment only one vertex x)
		 * g.removeVertex("a"); //Vertex x = opposite("a","first");//Should throw an
		 * error and print does not exist //Vertex x = opposite("c","third");//Should
		 * not get affected by the removing of c and print c Vertex x =
		 * opposite("c","second");//Should not get affected by the removing of c and
		 * print b
		 * 
		 * System.out.println(x._strUniqueID);
		 */

		/*
		 * //Test 4 testing edges method with and without removing vertix or edge
		 * Vector<Edge>v=new Vector<Edge>(); v=g.edges(); for(int i=0;i<v.size();i++) {
		 * System.out.println(v.get(i)._strUniqueID);
		 * System.out.println(v.get(i)._strData);
		 * System.out.println(v.get(i)._nEdgeCost);
		 * System.out.println("--------------"); }
		 * 
		 * System.out.
		 * println("-----After removing c the edge that is named third only 1 edge should be left ---------"
		 * ); g.removeVertex("c");
		 * 
		 * System.out.println("-----After removing edge named first ---------");
		 * g.removeEdge("first"); g.removeEdge("third"); v=g.edges(); for(int
		 * i=0;i<v.size();i++) { System.out.println(v.get(i)._strUniqueID);
		 * System.out.println(v.get(i)._strData);
		 * System.out.println(v.get(i)._nEdgeCost);
		 * System.out.println("--------------"); }
		 */

		/*
		 * //Test 5 printing Vertices Vector<Vertex> v=g.vertices(); for(int
		 * i=0;i<v.size();i++) { System.out.println(v.get(i)._strUniqueID);
		 * System.out.println(v.get(i)._strData); System.out.println(v.get(i)._nX);
		 * System.out.println(v.get(i)._nY); System.out.println("--------------"); }
		 * System.out.println("------ After removing Vertex b--------");
		 * g.removeVertex("b"); v=g.vertices(); for(int i=0;i<v.size();i++) {
		 * System.out.println(v.get(i)._strUniqueID);
		 * System.out.println(v.get(i)._strData); System.out.println(v.get(i)._nX);
		 * System.out.println(v.get(i)._nY); System.out.println("--------------"); }
		 */

		// Test 6 testing endVertices
		/*
		 * Vertex[] vert=endVertices("first"); //Vertex[]
		 * vert=endVertices("fifth");//should throw an error for(int
		 * i=0;i<vert.length;i++) { System.out.println(vert[i]._strUniqueID);
		 * System.out.println(vert[i]._strData); System.out.println(vert[i]._nX);
		 * System.out.println(vert[i]._nY); System.out.println("--------------"); }
		 */

		/*
		 * //TODO In test 7 it is throwing the wrong error in remove vertex however
		 * remove edge is working correctly //Test 7 removing vertex or edge that does
		 * not exist g.removeVertex("x"); //g.removeEdge("fifth");
		 * 
		 */
		/*
		 * //Test 8 Testing incident edges //changed in code because it was only working
		 * one sided but now it is correct //g.removeEdge("fourth");
		 * Vector<Edge>vctr=g.incidentEdges("a"); for(int i=0;i<vctr.size();i++) {
		 * System.out.println(vctr.get(i)._strUniqueID); }
		 */
	}
}
