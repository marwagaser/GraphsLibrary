import java.util.*; 
import java.lang.*; 
import java.io.*; 
  
class TesterGraph
{ 

  private int V;   // No. of vertices 
	    private LinkedList<Integer> tester[]; // testeracency List Represntation 
	  
	    // Constructor 
	    TesterGraph(int v) { 
	        V = v; 
	        tester = new LinkedList[v]; 
	        for(int i=0; i<v; ++i) 
	            tester[i] = new LinkedList(); 
	    } 
	  
	    // Function to add an edge into the graph 
	    boolean addEdge(String v1,String v2, TesterGraph g1) { 
	    	int v;
	    	int w;
	    	
	    	if (((int)v1.charAt(0))>=97&&((int)v1.charAt(0))<=122) {
	    		
	    		 v = ((int) v1.charAt(0))-97;
		    	 w = ((int) v2.charAt(0))-97;
	    	}
	    	
	    	else if (((int)v1.charAt(0))>=65&&((int)v1.charAt(0))<=90) {
	    		
	    		 v = ((int) v1.charAt(0))-65;
		    	 w = ((int) v2.charAt(0))-65;
	    	}
	    	else {
	    		 v = Integer.parseInt(v1);
	    		 w = Integer.parseInt(v2);
	    	}
	    	//System.out.println(w+" "+ v);
	    	//else if (((int)v1.charAt(0))>=97)
	        tester[v].add(w); 
	        tester[w].add(v); 
	        if (g1.isCyclic()) {
	        //	System.out.println("Cyclicccc");
	        	tester[v].removeLast();
		        tester[w].removeLast();
               return true;
	        }
	        else{
	        	
return false;
	        }
	    } 
	    
	 
	  
	    // A recursive function that uses visited[] and parent to detect 
	    // cycle in subgraph reachable from vertex v. 
	    Boolean isCyclicUtil(int v, Boolean visited[], int parent) 
	    { 
	        // Mark the current node as visited 
	        visited[v] = true; 
	        Integer i; 
	  
	        // Recur for all the vertices testeracent to this vertex 
	        Iterator<Integer> it = tester[v].iterator(); 
	        while (it.hasNext()) 
	        { 
	            i = it.next(); 
	  
	            // If an testeracent is not visited, then recur for that 
	            // testeracent 
	            if (!visited[i]) 
	            { 
	                if (isCyclicUtil(i, visited, v)) 
	                    return true; 
	            } 
	  
	            // If an testeracent is visited and not parent of current 
	            // vertex, then there is a cycle. 
	            else if (i != parent) 
	                return true; 
	        } 
	        return false; 
	    } 
	  
	    // Returns true if the graph contains a cycle, else false. 
	    Boolean isCyclic() 
	    { 
	        // Mark all the vertices as not visited and not part of 
	        // recursion stack 
	        Boolean visited[] = new Boolean[V]; 
	        for (int i = 0; i < V; i++) 
	            visited[i] = false; 
	  
	        // Call the recursive helper function to detect cycle in 
	        // different DFS trees 
	        for (int u = 0; u < V; u++) 
	            if (!visited[u]) // Don't recur for u if already visited 
	                if (isCyclicUtil(u, visited, -1)) 
	                    return true; 
	  
	        return false; 
	    } 
	  
	  
	    // Driver method to test above methods 
	   
	} 
	

