import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

public class Graph 
{

	// Keep a fast index to nodes in the map
	private Map<Integer, Vertex> vertexNames;

	/**
	 * Construct an empty Graph with a map. The map's key is the name of a vertex
	 * and the map's value is the vertex object.
	 */
	public Graph() 
	{
		vertexNames = new HashMap<>();
	}

	/**
	 * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
	 * with the same name are added.
	 * 
	 * @param v
	 *          (Vertex) vertex to be added to the graph
	 */
	public void addVertex(Vertex v) {
		if (vertexNames.containsKey(v.name))
			throw new IllegalArgumentException("Cannot create new vertex with existing name.");
		vertexNames.put(v.name, v);
	}

	/**
	 * Gets a collection of all the vertices in the graph
	 * 
	 * @return (Collection<Vertex>) collection of all the vertices in the graph
	 */
	public Collection<Vertex> getVertices() {
		return vertexNames.values();
	}

	/**
	 * Gets the vertex object with the given name
	 * 
	 * @param name
	 *          (String) name of the vertex object requested
	 * @return (Vertex) vertex object associated with the name
	 */
	public Vertex getVertex(String name) {
		return vertexNames.get(name); //using key to get vertex object
	}

	/**
	 * Adds a directed edge from vertex u to vertex v
	 * 
	 * @param nameU
	 *          (String) name of vertex u
	 * @param nameV
	 *          (String) name of vertex v
	 * @param cost
	 *          (double) cost of the edge between vertex u and v
	 */
	public void addEdge(int nameU, int nameV, Double cost) {
		//if the map doesn't contain the vertex then throw exceptions
		if (!vertexNames.containsKey(nameU))
			throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
		if (!vertexNames.containsKey(nameV))
			throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
		//use getVertex method to get the vertex object through the name
		Vertex sourceVertex = vertexNames.get(nameU);
		Vertex targetVertex = vertexNames.get(nameV);
		// create new objected directed from vertex to the other
		Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
		//adds it to source vertex's adjaceny listI suppose?
		sourceVertex.addEdge(newEdge);
	}

	/**
	 * Adds an undirected edge between vertex u and vertex v by adding a directed
	 * edge from u to v, then a directed edge from v to u
	 * 
	 * @param name
	 *          (String) name of vertex u
	 * @param name2
	 *          (String) name of vertex v
	 * @param cost
	 *          (double) cost of the edge between vertex u and v
	 */
	public void addUndirectedEdge(int name, int name2, double cost) { //undirected just means you add it both ways
		addEdge(name, name2, cost);
		addEdge(name2, name, cost);
	}


	/**
	 * Computes the euclidean distance between two points as described by their
	 * coordinates
	 * 
	 * @param ux
	 *          (double) x coordinate of point u
	 * @param uy
	 *          (double) y coordinate of point u
	 * @param vx
	 *          (double) x coordinate of point v
	 * @param vy
	 *          (double) y coordinate of point v
	 * @return (double) distance between the two points
	 */
	public double computeEuclideanDistance(double ux, double uy, double vx, double vy) { 
		return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2)); //uses x & y coordinates and computes distances through distance formula
	}

	/**
	 * Computes euclidean distance between two vertices as described by their
	 * coordinates
	 * 
	 * @param u
	 *          (Vertex) vertex u
	 * @param v
	 *          (Vertex) vertex v
	 * @return (double) distance between two vertices
	 */
	public double computeEuclideanDistance(Vertex u, Vertex v) { //returns the distance between two vertexes, just calls the other method
		return computeEuclideanDistance(u.x, u.y, v.x, v.y);
	}

	/**
	 * Calculates the euclidean distance for all edges in the map using the
	 * computeEuclideanCost method.
	 */
	public void computeAllEuclideanDistances() { //what I coded before, computes all distances for all edges between vertices
		for (Vertex u : getVertices()) //traverse through vertex objects 
			for (Edge uv : u.adjacentEdges) { //traverse through adjacency list of vertices
				Vertex v = uv.target; //find target of edge (where its pointing to, its source is u)
				uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y); //compute distance and set it as edge length
			}
	}



	// STUDENT CODE STARTS HERE

	public void generateRandomVertices(int n) {
		//  private Map<Integer, Vertex> vertexNames;
		vertexNames = new HashMap<>(); // reset the vertex hashmap

		//generate n ranadom vertices
		for(int i = 0; i < n; i++) 
		{
			//so name of vertex would be 0, 1, 2, 3... n - 1
			//new vertex with index as name, and generated random #s from 0-100 (inclusive) as x and y coordinates
			Vertex v = new Vertex(i, (int)(Math.random()*101), (int)(Math.random()*101)); 

			//add it to the map
			vertexNames.put(i, v); //i is the name or key and v is the vertex object or the value

		}

		//ok so now going through, put undirected edges iterating through at each index and beyond
		for(int i = 0; i < n; i++)
		{
			for(int j = i+1; j < n; j++)
			{
				Vertex mainVertex = vertexNames.get(i);
				Vertex otherVertex = vertexNames.get(j);
				addUndirectedEdge(i, j, computeEuclideanDistance(mainVertex.x, mainVertex.y, otherVertex.x, otherVertex.y));
			}
		}

		//do we set the vertex adjacaceny lists? I don't understand how this method can happen?
		//I think we set every single one of the vertex's adjacency edge to be every other vertex but itself

		//is this even needed? shoudln't it set automatically?
		computeAllEuclideanDistances(); // compute distances and set all edges to correct distances
		//does it twice but we can just leave it 
	}

	public List<Edge> nearestNeighborTsp() 
	{
		//pick a random starting city by randomizing key and getting vertex of it in the map
		//since indexes would be 0, 1, 2, 3, 4, .... n-1, then math.random * n (which is the map size)
		Vertex start = vertexNames.get((int)(Math.random()*vertexNames.size()));
		System.out.println("random city is " + start.name);

		//initialize list of edges representing shortest paths of city to each other
		List<Edge> shortest = new LinkedList<Edge>();

		//list of viisted cities to keep track of in order to not going to them again
		List<Vertex> visited = new LinkedList<Vertex>();
		visited.add(start);

		System.out.println("Adding to not visited: ");
		//list of all vertices that have not been visited yet
		List<Vertex> notVisited = new LinkedList<Vertex>();
		//putting all vertices inside (based on index)
		for(int i = 0; i < vertexNames.size(); i++)
		{
			notVisited.add(vertexNames.get(i)); //now notVisited list contains all the vertices
			System.out.println("added to notVistied: " + i);
		}

		notVisited.remove(start.name); //remove the starting one from the list because it's visited (index same as name in list)

		//after adding all cities the size of shortest would be the map size so you can stop searching
		while(!notVisited.isEmpty())
		{
			//iterate through the adjacency list of start and find shortest path

			//set starting min index to first one that is not in visited

			//out of all that are not visited, set current min to the first ones
			//then iterate through and find the actual min, then remove it

			//temp min is beginning from start to first vertex in not visited list 
			Edge min = edgeDistanceBetweenTwoVertices(start, notVisited.get(0));
			int minIndex = 0; //index of the vertex name of the target of current vertex's edge

			for(int i = 1; i < notVisited.size(); i++)
			{
				Edge e = edgeDistanceBetweenTwoVertices(start, notVisited.get(i));
				if(e.distance < min.distance) 
				{
					min = e;
					minIndex = i;
				}
			}

				//should end with the min one, add edge to list
				shortest.add(min);
				
				System.out.println("added min edge from " + min.source + " to " + min.target + " with length of " + min.distance);

				System.out.println("shortest size is: " + shortest.size());
				//don't know if visited even maters honestly, but for now add to visited list
				visited.add(min.target);
				
				System.out.println("added " + min.target + " to visited vertices");

				start = min.target;
				
				System.out.println("The new start is " + start);

				notVisited.remove(min.target); //remove it form the not visited list since now you have visited it

			



			//			for(Edge e: visited.get(visited.size()-1).adjacentEdges)
			//			{
			//				if (e.target == visited.get(0)) //if target of last one edge points back to first one
			//				{
			//					//then this is the correct edge, add it to shortest list
			//					shortest.add(e);
			//				}
			//			}

		}
		
		//add the last edge that points to first one back
		//the last 1 of visited is the last city
		//find edge that goes back to the first one
		Edge endToStart = edgeDistanceBetweenTwoVertices(visited.get(visited.size()-1),visited.get(0));
		shortest.add(endToStart);
		
		//calcualte total distance for nearest neight
		int sum = 0;
		for(int i = 0; i < shortest.size(); i++)
		{
			sum += shortest.get(i).distance;
		}
		
		System.out.println("NEAREST NEIGHBOR SUM IS: "+ sum);
		
		return shortest; //return shortest which is a list of edges that points from the source to target
	}

	//returns the edge between two vertices
	public Edge edgeDistanceBetweenTwoVertices(Vertex source, Vertex target)
	{
		for(Edge e: source.adjacentEdges)
		{
			if(e.target == target)
			{
				return e;
			}
		}

		return null;
	}

	public List<Edge> bruteForceTsp() 
	{
		//Enumerate all possible permutations of the integers 0 to n-1
		ArrayList<Integer> finalPermutations = new ArrayList<Integer>();

		int[] indexesArr = new int[vertexNames.size()];

		//first make an array of indexes from 0 to n-1
		for(int i = 0; i < vertexNames.size(); i++)
		{
			indexesArr[i] = i;
		}

		//arraylist of the permuted indexes stored as strings
		ArrayList<String>  strOfPermutations = new ArrayList<String>();

		//function to permute the indexes giving them the start, array of indexes, returns arraylist of permuted indexes in string format
		permute(0, indexesArr, strOfPermutations);    
		
//		System.out.println("String of permutations: ");
//		for(int i = 0; i < strOfPermutations.size(); i++)
//		{
//			System.out.println(strOfPermutations.get(i));
//
//		}

		//List to store all sum values of each path/permutation
		LinkedList<Integer> sumOfPermutationPaths = new LinkedList<Integer>();

		//for each string permutation, split each number into a different index in an arraylist
		for(String x: strOfPermutations)
		{			
			//get the edge list of each permutation and check its total distance
			List<Edge> edgeListMap = getEdgeListGivenPermutationString(x);
			//System.out.println("edgeList for strings: " + edgeListMap);
			//iterate through edge list and get sum
			//set initial sum to 0
			int sum = 0;
			
			for(int i = 0; i < edgeListMap.size(); i++)
			{
				sum += edgeListMap.get(i).distance;
			}
			
			sumOfPermutationPaths.add(sum);
			
		}
		
//		for(int i = 0; i < sumOfPermutationPaths.size(); i++)
//		{
//			System.out.println("Sum for iteration " + i + " of " + strOfPermutations.get(i) + " order is " + sumOfPermutationPaths.get(i));
//
//		}

		//index of the min	-- set this as 0 min for now
		int indexOfMinPath = 0;
		int min = sumOfPermutationPaths.get(0);

		//iterate through the sum of permutation path and find the smallest
		for(int i = 1; i < sumOfPermutationPaths.size(); i++)
		{
			if(sumOfPermutationPaths.get(i) < min)
			{
				indexOfMinPath = i;
				min = sumOfPermutationPaths.get(i);
			}

		}
		
//		System.out.println("The min index is " + indexOfMinPath);
//		
//		System.out.println("The min order is " + strOfPermutations.get(indexOfMinPath));
//		
//		System.out.println("BRUTE FORCE TSP MIN SUM IS " + sumOfPermutationPaths.get(indexOfMinPath));


		List<Edge> finalMinPath = getEdgeListGivenPermutationString(strOfPermutations.get(indexOfMinPath));

		

		return finalMinPath; 
	}

	//get permutation of edges in a list given string representation
	public List<Edge> getEdgeListGivenPermutationString(String s)
	{
		//for each string permutation, split each number into a different index in an arraylist

		ArrayList<Integer> x = (ArrayList<Integer>) Arrays.stream(s.split("\\s"))
				.map(Integer::parseInt)
				.collect(Collectors.toList());
		
		//System.out.println("x is"  + x);

		//then create a list that contains the edges of each
		List<Edge> permutedEdges = new LinkedList<Edge>();
		
		for(int k = 0; k < x.size(); k++)
		{
			//if its the last index, get edge from last to start and add
			if(k == x.size()-1)
			{
				permutedEdges.add(edgeDistanceBetweenTwoVertices(vertexNames.get(x.get(k)), vertexNames.get(x.get(0))));
			}
			else
			{
				permutedEdges.add(edgeDistanceBetweenTwoVertices(vertexNames.get(x.get(k)), vertexNames.get(x.get(k+1))));

			}


		}
		
		return permutedEdges;
	}

		//Enumerate all possible permutations of the integers 0 to n-1
		public static void permute(int start, int[] input, ArrayList<String> list)
		{

			if (start == input.length) {
				String hello = "";
				//System.out.println(input);
				for(int x: input)
				{

					//System.out.print(x + " ");
					hello += x + " ";

					//return x;
				}
				//System.out.print("This is hello: " + hello);
				//System.out.println("");
				list.add(hello);
				//System.out.println(count);
				//return"";
			}
			for (int i = start; i < input.length; i++) {
				// swapping
				int temp = input[i];
				input[i] = input[start];
				input[start] = temp;
				// swap(input[i], input[start]);

				permute(start + 1, input, list);
				// swap(input[i],input[start]);

				int temp2 = input[i];
				input[i] = input[start];
				input[start] = temp2;
			}
		}

		// STUDENT CODE ENDS HERE



		/**
		 * Prints out the adjacency list of the graph for debugging
		 */
		public void printAdjacencyList() {
			for (int u : vertexNames.keySet()) {
				StringBuilder sb = new StringBuilder();
				sb.append(u);
				sb.append(" -> [ ");
				for (Edge e : vertexNames.get(u).adjacentEdges) {
					sb.append(e.target.name);
					sb.append("(");
					sb.append(e.distance);
					sb.append(") ");
				}
				sb.append("]");
				System.out.println(sb.toString());
			}
		}
	}
