package program314;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;


public class Deploy
{	
	private static final int Infinity = Integer.MAX_VALUE/2;
	
	private static List<String> result;
	
	private static List<Integer> minSource;
	
	private static int minCost = Infinity;
	
	private static int[][] minFlow;
	
    private static int numOfConsumers;

	private static List<Consumer> consumers;
	
	/**
	 * 由consumer.link查找consumer.id
	 */
	private static Map<Integer, Integer> getConsumerId;

	private static int serverCost;

	private static int numOfVertice;
	
	private static int diNumOfVertice;
	
	private static int superSource;
	
	private static int superSink;
	
	/**
     * 你需要完成的入口
     * <功能详细描述>
     * @param graphContent 用例信息文件
     * @return [参数说明] 输出结果信息
     * @see [类、类#方法、类#成员]
     */
    public static String[] deployServer(String[] graphContent) {
    	//原始网络
    	GraphO graphO = createGraphO(graphContent);
    	//原始网络的有向边
    	List<List<GraphO.Edge>> edges = graphO.edges;
    	
    	//消费节点i的id和顶点j的下标的有效距离
    	int[][] distance = new int[numOfConsumers][numOfVertice];
		shortestDistance(edges, distance);
		
		
		//计时器
   		Time time = new Time();
   		int maxTime = 60000;
		for (int i = 0; i < 1000; i++) {
			if (time.getTimeDelay() > maxTime) {
				break;
			}
			
			//获取源点的组合
	    	List<List<Integer>> situation = new ArrayList<>(numOfConsumers);
			getSourceSituation(distance, situation);
			
			//计算该源点组合下的最小费用流
	    	for (List<Integer> sources : situation) {
	    		getMinFlow(graphO, sources);
	    	}
	    }
		
		//计算极端情况,所有消费节点连接服务器
		List<Integer> extreme = new ArrayList<>();
		for (Consumer consumer : consumers) {
			extreme.add(consumer.link);
		}
		getMinFlow(graphO, extreme);
    	
		//caseEx.txt的最优解?
//    	minSource = new ArrayList<>(Arrays.asList(new Integer[]{0, 3, 22}));
//    	minFlow = new int[][]{{0,20,16,0,0,0,13,10,36,13,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,7,0,11,2,0,0,0,0,0,0,0,0},{0,2,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0,0},{0,11,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,12,0,0,0},{0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,11,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    	
	
		
    	getResult(graphO);
       	
    	return result.toArray(new String[result.size()]);
    }
	/**
	 * 通过最小费用流矩阵求得路径
	 * @param graphO
	 * @param edges
	 */
	private static void getResult(GraphO graphO) {
		
		List<List<GraphO.Edge>> edges = graphO.edges; 
		
		//修改原始网络和minFlow，增加虚拟源汇点和对应的边
    	//虚拟源点的邻接表
    	List<GraphO.Edge> sourceEdges = new ArrayList<>(minSource.size());
    	for (int source : minSource) {
			sourceEdges.add(new GraphO.Edge(numOfVertice, source, Infinity, 0));
			graphO.numOfEdges++;
		}
    	edges.add(sourceEdges);
    	//虚拟汇点的邻接表
    	for (Consumer consumer : consumers) {
			edges.get(consumer.link).add(new GraphO.Edge(consumer.link, numOfVertice + 1, consumer.demand, 0));
			graphO.numOfEdges++;
		}
    	edges.add(new ArrayList<GraphO.Edge>());
    	//最小费用流矩阵
    	int[][] temp = minFlow;
    	minFlow = new int[numOfVertice + 2][numOfVertice + 2];
    	for (int i = 0; i < numOfVertice; i++) {
			for (int j =0; j < numOfVertice; j++) {
				minFlow[i][j] = temp[i][j];
			}
		}
    	for (int source : minSource) {
			minFlow[numOfVertice][source] = Infinity;
		}
    	for (Consumer consumer : consumers) {
			minFlow[consumer.link][numOfVertice + 1] = consumer.demand;
		}
    	
    	result = new ArrayList<>();

    	while (true) {
			if (!dfs(graphO, new boolean[numOfVertice + 2], numOfVertice, Infinity, new ArrayList<Integer>())) {
				break;
			}
		}
    	int count = result.size();
    	
    	result.addAll(0, Arrays.asList(new String[]{Integer.toString(count), ""}));
	}
    /**
     * @param graphO
     * @param visited
     * @param cur
     * @param flow 源点到cur的路径流量
     * @param path
     * @return 是否找到虚拟源点到虚拟汇点的路径
     */
    private static boolean dfs(GraphO graphO, boolean[] visited, int cur, int flow, List<Integer> path) {
    	if (cur == numOfVertice + 1) {
    		path.add(cur);
			StringBuilder builder = new StringBuilder();
			for (int i = 1; i < path.size() - 1; i++) {
				builder.append(path.get(i) + " ");
			}
			int consumerId = getConsumerId.get(path.get(path.size() - 2));
			builder.append(consumerId + " ");
			builder.append(flow);
			
//			logger.println(builder.toString());
//			logger.flush();
			
			result.add(builder.toString());
			
			for (int i = 0; i < path.size() - 1; i++) {
				minFlow[path.get(i)][path.get(i+1)] -= flow;
			}
			return true;
		}
    	
    	path.add(cur);
    	visited[cur] = true;
    	
    	for (GraphO.Edge edge : graphO.edges.get(cur)) {
			int next = edge.end;
			if (!visited[next] && minFlow[cur][next] > 0) {
				if (dfs(graphO, visited, next, Math.min(flow, minFlow[cur][next]), path)) {
					return true;
				}
			}
		}
    	path.remove(path.size() - 1);
    	visited[cur] = false;
    	
    	return false;
    }
	/**
	 * @param graphO
	 * @param result
	 * @param sources
	 */
	private static void getMinFlow(GraphO graphO, List<Integer> sources) {
		GraphDi graphDi = new GraphDi(graphO, sources);
		
		int[][] flowDi = new int[diNumOfVertice][diNumOfVertice];
		int[][] remain = new int[diNumOfVertice][diNumOfVertice];
		
		//初始化剩余网络，添加反向弧
		for (List<Arc> links : graphDi.arcs) {
			for (Arc arc : links) {
				if (arc.isback) {
					continue;
				}
				remain[arc.start][arc.end] = arc.capacity;
				Arc backArc = new Arc(arc.end, arc.start, arc.capacity, -arc.cost);
				backArc.isback = true;
				graphDi.arcs.get(backArc.start).add(backArc);
			}
		}
		
		while (true) {
			//求superSource到superSink的最短路径
			Path[] paths = new Path[diNumOfVertice];
			for (Arc arc : graphDi.arcs.get(superSource)) {
				paths[arc.end] = new Path(arc.end, superSource, arc.cost);
			}
			for (int vertex = 0; vertex < paths.length; vertex++) {
				if (vertex == superSource) {
					//pre=-2:起点
					paths[vertex] = new Path(vertex, -2, 0);
				} else if (paths[vertex] == null) {
					paths[vertex] = new Path(vertex, -1, Infinity);
				}
			}
			
			for (int vertex = 0; vertex < paths.length; vertex++) {
				for (List<Arc> links : graphDi.arcs) {
					for (Arc arc : links) {
						int start = arc.start;
						int end = arc.end;
						int cost = arc.cost;
						if (paths[start].pre != -1 && remain[start][end] != 0) {
							if (paths[start].distance + cost < paths[end].distance) {
								paths[end].distance = paths[start].distance + cost;
								paths[end].pre = start;
							}
						}
					}
				}
			}
			
			//从superSink到superSource的最短路径
			List<Integer> shortestPath = new ArrayList<>();
			int cur = superSink;
			while (cur >= 0) {
				shortestPath.add(cur);
				cur = paths[cur].pre;
			}
			
			//无最短路径
			if (cur == -1) {
				break;
			}
			
			//最短路径的流量
			int minFlow = Infinity;
			for (int i = shortestPath.size() - 1; i >= 1; i--) {
				int start = shortestPath.get(i);
				int end = shortestPath.get(i-1);
				minFlow = Math.min(minFlow, remain[start][end]);
			}
			
//				输出有向图最短路径
//				writer.println("flow:"+minFlow);
//				writer.print("shortestPath:");
//				for (int i = shortestPath.size() - 1; i >= 0; i--) {
//					writer.print(shortestPath.get(i)+" ");
//				}
//				writer.println();
			
			//求剩余网络
			for (int i = shortestPath.size() - 1; i >= 1; i--) {
				int start = shortestPath.get(i);
				int end = shortestPath.get(i-1);
				Arc section = null;
				for (Arc arc : graphDi.arcs.get(start)) {
					if (arc.start == start) {
						section = arc;
						break;
					}
				}
				//section是否为反向弧
				if (!section.isback) {
					flowDi[start][end] += minFlow;
					remain[start][end] -= minFlow;
					remain[end][start] = flowDi[start][end];
				} else {
					flowDi[start][end] -= minFlow;
					remain[start][end] += minFlow;
					remain[end][start] = flowDi[start][end];
				}
			}
		}
		
		//判断是否满足消费节点的需求
		boolean enough = true;
		for (Consumer consumer : consumers) {
			int sink = consumer.link + numOfVertice;
			if (flowDi[sink][superSink] != consumer.demand) {
				enough = false;
			}
//				输出每个消费节点的需求和实际流量
//				writer.println("costumerout "+sink+" need "+consumer.demand+
//						"have " + flowDi[sink][superSink]);
		}
//		writer.println("isEnough:"+enough);
		if (enough) {
			//转化为无向图中的flow
			int[][] flow = new int[numOfVertice][numOfVertice];
			for (int i = numOfVertice; i < 2 * numOfVertice; i++) {//出节点
				for (int j = 0; j < numOfVertice; j++) {//入节点
					flow[i - numOfVertice][j] = flowDi[i][j];
				}
			}
			int[][] cost = new int[numOfVertice][numOfVertice];
			for (List<GraphO.Edge> links : graphO.edges) {
				for (GraphO.Edge edge : links) {
					cost[edge.start][edge.end] = edge.cost;
				}
			}
			
			int sumCost = sources.size() * serverCost;
			
			for (int i = 0; i < numOfVertice; i++) {
				for (int j = 0; j < i; j++) {
					sumCost += (Math.max(flow[i][j], flow[j][i])) * cost[i][j];
				}
			}
			
			if (sumCost < minCost) {
				minCost = sumCost;
				minFlow = flow;
				minSource = sources;
			}
		}
		
	}
    
    /**
	 * k-means聚类寻找源点
	 * @param distance
	 * @param situation
	 */
    
	private static void getSourceSituation(int[][] distance, List<List<Integer>> situation) {
		
    	//初始num个源点时的情况
    	for (int num = 1; num < numOfConsumers; num++) {
    		//源点的位置
			List<Integer> position = new ArrayList<>(num);
			while (position.size() < num) {
				int temp = (int) (Math.random() * numOfVertice);
				if (!position.contains(temp)) {
					position.add(temp);
				}
			}
			//K=源点下标 V=消费节点id的集合
			Map<Integer, List<Integer>> map = new HashMap<>(num);
			for (int consumerId = 0; consumerId < numOfConsumers; consumerId++) {
				int minSource = -1;
				int minDistance = Infinity;
				for (int source : position) {
					if (distance[consumerId][source] < minDistance) {
						minSource = source;
						minDistance = distance[consumerId][minSource];
					}
				}
				if (!map.containsKey(minSource)) {
					map.put(minSource, new ArrayList<Integer>());
				}
				List<Integer> consumerCollection = map.get(minSource);
				consumerCollection.add(consumerId);
			}
			boolean move = true;
			while (move) {
				move = false;
				Map<Integer, List<Integer>> newMap = new HashMap<>(num);
				for (int source : map.keySet()) {
					List<Integer> consumerCollection = map.get(source);
					//新source对应的消费节点集合的中心
					int minCenter = -1;
					int minSum = Infinity;
					for (int vertex = 0; vertex < numOfVertice; vertex++) {
						int sum = 0;
						for (int consumerId : consumerCollection) {
							sum += distance[consumerId][vertex];
						}
						if (sum < minSum) {
							minCenter = vertex;
							minSum = sum;
						}
					}
					if (source != minCenter) {
						move = true;
					} 
					if (newMap.containsKey(minCenter)) {
						newMap.get(minCenter).addAll(consumerCollection);
					} else {
						newMap.put(minCenter, consumerCollection);
					}
				}
				map = newMap;
			}
			
			situation.add(new ArrayList<>(map.keySet()));
		}
//    	writer.println("situation[ num:"+situation.size()+" content:"+situation.toString()+"]");
	}

	/**
	 * 求原始网络中汇点到顶点的最短距离
	 * @param edges
	 * @param distance
	 */
	private static void shortestDistance(List<List<GraphO.Edge>> edges, int[][] distance) {
		for (int situation = 0; situation < numOfConsumers; situation++) {
    		Consumer consumer = consumers.get(situation);
			int sink = consumer.link;
			
			Path[] paths = dijstra(edges, sink);
			
			for (int i = 0; i < numOfVertice; i++) {
				distance[situation][i] = paths[i].distance;
			}
		}
	}

	/**
	 * @param edges
	 * @param v0 最短路径的起点
	 * @return
	 */
	private static Path[] dijstra(List<List<GraphO.Edge>> edges, int v0) {
		Path[] paths = new Path[numOfVertice];
		for (int vertex = 0; vertex < numOfVertice; vertex++) {
			if (vertex == v0) {
				paths[vertex] = new Path(vertex, -2, 0);
			} else {
				int dis = getWeight(v0, vertex, edges);
				if (dis == Infinity) {
					paths[vertex] = new Path(vertex, -1, dis);
				} else {
					paths[vertex] = new Path(vertex, v0, dis);
				}
			}
		}
		
		//存储没有加入最短路径的节点
		PriorityQueue<Path> heap = new PriorityQueue<>(numOfVertice, new Comparator<Path>() {

			@Override
			public int compare(Path o1, Path o2) {
				return o1.distance-o2.distance;
			}
		});
		for (int i = 0; i < numOfVertice; i++) {
			if (i != v0) {
				heap.add(paths[i]);
			}
		}
		
		while (!heap.isEmpty()) {
			Path cur = heap.poll();
			int start = cur.id;
			for (GraphO.Edge edge : edges.get(start)) {
				int next = edge.end;
				if (paths[next].pre == -1 && 
						cur.distance + edge.cost < paths[next].distance) {
					heap.remove(paths[next]);
					paths[next].pre = start;
					paths[next].distance = cur.distance + edge.cost;
					heap.add(paths[next]);
				}
			}
		}
		return paths;
	}

	private static int getWeight(int start, int end, List<List<GraphO.Edge>> edges) {
    	List<GraphO.Edge> links = edges.get(start);
    	if (links.isEmpty()) {
			return Infinity;
		} else {
			for (GraphO.Edge edge : links) {
				if (edge.end == end) {
					return edge.cost;
				}
			}
			return Infinity;
		}
    }
    private static GraphO createGraphO(String[] graphContent) {
		String[] tokens = graphContent[0].split(" ");
    	int numOfVertice = Integer.parseInt(tokens[0]);
    	int numOfEdges = Integer.parseInt(tokens[1]);
    	int numOfConsumers = Integer.parseInt(tokens[2]);
    	int serverCost = Integer.parseInt(graphContent[2]);
    	
    	Deploy.numOfVertice = numOfVertice;
		Deploy.diNumOfVertice = numOfVertice * 2 + 2;
//		Deploy.pNumOfVertice = diNumOfVertice + 2;
    	Deploy.numOfConsumers = numOfConsumers;
    	Deploy.superSource = 2 * numOfVertice;
		Deploy.superSink = superSource + 1;
//		Deploy.plusSource = superSink + 1;
//		Deploy.plusSink = plusSource + 1;
		
    	Deploy.serverCost = serverCost;
    	Deploy.numOfConsumers = numOfConsumers;
    	
    	GraphO graphO = new GraphO(numOfEdges);
    	List<List<GraphO.Edge>> edges = new ArrayList<>(numOfVertice);
    	for (int i = 0; i < numOfVertice; i++) {
			edges.add(new ArrayList<GraphO.Edge>());
		}
    	graphO.edges = edges;
    	
    	int index;
    	for (index = 4; index < graphContent.length; index++) {
    		String line = graphContent[index];
    		if (line.length()==0) {
				break;
			}
			String[] edge = line.split(" ");
			int v1 = Integer.parseInt(edge[0]);
			int v2 = Integer.parseInt(edge[1]);
			int capacity = Integer.parseInt(edge[2]);
			int cost = Integer.parseInt(edge[3]);
			edges.get(v1).add(new GraphO.Edge(v1, v2, capacity, cost));
			edges.get(v2).add(new GraphO.Edge(v2, v1, capacity, cost));
		}
    	
    	Deploy.consumers = new ArrayList<>(numOfConsumers);
    	Deploy.getConsumerId = new HashMap<>(numOfConsumers);
    	
    	for (index++; index < graphContent.length; index++) {
			String[] consumer = graphContent[index].split(" ");
			int id = Integer.parseInt(consumer[0]);
			int link = Integer.parseInt(consumer[1]);
			int demand = Integer.parseInt(consumer[2]);
			consumers.add(new Consumer(id, link, demand));
			getConsumerId.put(link, id);
		}
    	
    	return graphO;
	}
    
    /**
	 * 原始无向图
	 * @author hjg
	 *
	 */
	private static class GraphO {
		
		private int numOfEdges;
		
		private List<List<Edge>> edges;
		
		
		public GraphO(int numOfEdges) {
			this.numOfEdges = numOfEdges;
		}
		
		private static class Edge {
			
			private int start;
			
			private int end;
			
			private int capacity;
			
			private int cost;
			
			public Edge(int start, int end, int capacity, int cost) {
				this.start = start;
				this.end = end;
				this.capacity = capacity;
				this.cost = cost;
			}
		}
	}

	private static class GraphDi {
		
		private List<List<Arc>> arcs;
		
		public GraphDi(GraphO graphO, List<Integer> sources) {
			List<List<Arc>> arcs = new ArrayList<>(diNumOfVertice);
			this.arcs = arcs;
			
			for (int i = 0; i < diNumOfVertice; i++) {
				arcs.add(new ArrayList<Arc>(5));
			}
			for (int vertex = 0; vertex < numOfVertice; vertex++) {
				arcs.get(vertex).add(new Arc(vertex, vertex + numOfVertice, Infinity, 0));
			}
			for (int vertex = 0; vertex < numOfVertice; vertex++) {
				for (GraphO.Edge edge : graphO.edges.get(vertex)) {
					int next = edge.end;
					if (vertex < next) {
						int in1 = vertex;
						int out1 = vertex + numOfVertice;
						int in2 = next;
						int out2 = next + numOfVertice;
						arcs.get(out1).add(new Arc(out1, in2, edge.capacity, edge.cost));
						arcs.get(out2).add(new Arc(out2, in1, edge.capacity, edge.cost));
					}
				}
			}
			List<Arc> superSourceLinks = arcs.get(superSource);
			for (int source : sources) {
				int inSource = source;
				superSourceLinks.add(new Arc(superSource, inSource, Infinity, 0));
			}
			for (Consumer consumer : consumers) {
				int outSink = consumer.link + numOfVertice;
				arcs.get(outSink).add(new Arc(outSink, superSink, consumer.demand, 0));
			}
		}
	}

	private static class Consumer {
		
		private int id;
		
		private int link;
		
		private int demand;
		
		public Consumer(int id, int link, int demand) {
			this.id = id;
			this.link = link;
			this.demand = demand;
		}
		
		@Override
		public String toString() {
			return "[id=" + id + " link=" + link + " demand=" + demand;
		}
	}

	private static class Path {
		
		private int id;
		
		private int pre;
		
		private int distance;
		
		public Path(int id, int pre, int distance) {
			this.id = id;
			this.pre = pre;
			this.distance = distance;
		}
	}

	private static class Arc {
		
		/**
		 * 是否反向弧
		 */
		private boolean isback;
		
		private int start;
		
		private int end;
		
		private int capacity;
		
		private int cost;
		
		public Arc(int start, int end, int capacity, int cost) {
			this.start = start;
			this.end = end;
			this.capacity = capacity;
			this.cost = cost;
			isback = false;
		}
		
		@Override
		public String toString() {
			return "["+start+"->"+end+": capacity "+capacity+" cost "+cost+" isBack "+isback+"]";
		}
	}

	public static class Time
	{
	    private static final long start = System.currentTimeMillis();
	
	    private long current = 0;
	
	    public Time()
	    {
	    }
	
	    public long getTimeDelay()
	    {
	        current = System.currentTimeMillis();
	        return current - start;
	    }
	
	    public long getStart()
	    {
	        return start;
	    }
	}
}
