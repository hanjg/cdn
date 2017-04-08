package program316;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;

import com.filetool.util.LogUtil;


public class Deploy
{	
	private static PrintWriter logger;
	
	static {
		try {
			logger = new PrintWriter(new File("log.txt"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private static final int Infinity = Integer.MAX_VALUE/2;
	
	private static int defaultSize;
	
	/**
	 * 最终输出的路径数量和路径
	 */
	private static List<String> result;
	
	private static List<Integer> minSource;
	
	private static int minCost = Infinity;
	
	private static int[][] minFlow;
	
	/**
	 * 有向图中除了超级源汇点之外的弧的费用
	 */
	private static int[][] cost;
	
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
	 * 有向图拆点拆边，添加超级汇点和对应边的弧
	 */
	private static List<List<Arc>> arcsDi;
	

	//计时器
	private static Time timer;
	/**
     * 你需要完成的入口
     * <功能详细描述>
     * @param graphContent 用例信息文件
     * @return [参数说明] 输出结果信息
     * @see [类、类#方法、类#成员]
     */
    public static String[] deployServer(String[] graphContent) {
    	timer = new Time();
   		int maxTime = 80000;
   		
   		LogUtil.printLog("startCore ");
   		
    	//原始网络
    	GraphO graphO = createGraphO(graphContent);
    	
    	//消费节点(可以看做和汇点重合)和顶点的有效距离(i为消费节点id，j为顶点下标)
    	int[][] distance = new int[numOfConsumers][numOfVertice];
		shortestDistance(graphO, distance);
		
		//输出网络的基本信息
		logger.println("numOfVertice:" + numOfVertice +
				"\nnumOfEdges:" + graphO.numOfEdges + 
				"\nnumOfConsumers:" + numOfConsumers);
				
		LogUtil.printLog("start get minFlow ");
		
		for (int i = 0; i < 10000; i++) {
			if (timer.getTimeDelay() > maxTime) {
				break;
			}
			//获取源点的组合
			List<List<Integer>> situation = getSourceSituation(distance);
			
			//计算该源点组合下的最小费用流
	    	for (List<Integer> sources : situation) {
	    		getMinFlow(graphO, sources);
	    	}
	    }
		
		//计算极端情况的费用：所有消费节点连接服务器
		List<Integer> extreme = new ArrayList<>(numOfConsumers);
		for (Consumer consumer : consumers) {
			extreme.add(consumer.link);
		}
		getMinFlow(graphO, extreme);
    			
		//caseEx.txt的最优解?
//    	minSource = new ArrayList<>(Arrays.asList(new Integer[]{0, 3, 22}));
//    	minFlow = new int[][]{{0,20,16,0,0,0,13,10,36,13,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,7,0,11,2,0,0,0,0,0,0,0,0},{0,2,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0,0},{0,11,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,12,0,0,0},{0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,11,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    	
		logger.println("minCost:" + minCost);
		logger.println("minSource:" + minSource);
		
		LogUtil.printLog("start getResult ");
		
    	getResult(graphO);
       	
    	logger.print("coreTime:" + timer.getTimeDelay() + " ms");
		
		logger.flush();
    	
    	return result.toArray(new String[result.size()]);
    }
    
	/**
	 * 通过最小费用流矩阵求得路径
	 * @param graphO
	 * @param edges
	 */
	private static void getResult(GraphO graphO) {
		List<List<GraphO.Edge>> edges = graphO.edges; 
		
		//修改原始网络，增加虚拟源汇点和对应的边
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
    	
		//有向图中的minFlow转化为无向图中的flow
		int[][] flow = new int[numOfVertice][numOfVertice];
		for (int i = numOfVertice; i < 2 * numOfVertice; i++) {//出节点
			for (int j = 0; j < numOfVertice; j++) {//入节点
				flow[i - numOfVertice][j] = minFlow[i][j];
			}
		}
    	int[][] temp = flow;
    	
    	//原始网络中添加超级源汇点，得到新的流网络
    	//复制原始网络的边
    	minFlow = new int[numOfVertice + 2][numOfVertice + 2];
    	for (int i = 0; i < numOfVertice; i++) {
			for (int j =0; j < numOfVertice; j++) {
				minFlow[i][j] = temp[i][j];
			}
		}
    	
    	//在流网络中增加超级源汇点对应的边
    	for (int source : minSource) {
			minFlow[numOfVertice][source] = Infinity;
		}
    	for (Consumer consumer : consumers) {
			minFlow[consumer.link][numOfVertice + 1] = consumer.demand;
		}
    	
    	result = new ArrayList<>(defaultSize);

    	while (true) {
			if (!dfs(graphO, new boolean[numOfVertice + 2], numOfVertice, Infinity, new ArrayList<Integer>(defaultSize))) {
				break;
			}
		}
    	int count = result.size();
    	
    	result.addAll(0, Arrays.asList(new String[]{Integer.toString(count),""}));
	}
	
	@SuppressWarnings("unused")
	private static void printMinFlow() {
		logger.print("minFlow:\n{");
		for (int i = 0; i < numOfVertice; i++) {
			if (i != 0) {
				logger.print(",");
			}
			logger.print("{");
			for (int j =0; j < numOfVertice; j++) {
				if (j != 0) {
					logger.print(",");
				}
				logger.print(minFlow[i][j]);
			}
			logger.print("}");
		}
		logger.print("}");
		logger.flush();
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
		//根据初始无向网络和源点位置，拆点，拆边，增加超级源汇点和对应弧，新建有向图
		GraphDi graphDi = new GraphDi(graphO, sources);
		
		//有向图的流网络
		int[][] flowDi = new int[diNumOfVertice][diNumOfVertice];
		//有向图的剩余网络
		int[][] remain = new int[diNumOfVertice][diNumOfVertice];
		
		//初始化流网络
		for (List<Arc> links : graphDi.arcs) {
			for (Arc arc : links) {
				if (!arc.isback) {
					remain[arc.start][arc.end] = arc.capacity;
				}
			}
		}
		
		while (true) {
			//从superSink到superSource的最短路径
			List<Integer> shortestPath = new ArrayList<>(defaultSize);
			
			//SPFA计算最短路径
			Path[] paths = new Path[diNumOfVertice];
			for (int vertex = 0; vertex < diNumOfVertice; vertex++) {
				if (vertex == superSource) {
					paths[vertex] = new Path(vertex, -2, 0);
				} else {
					paths[vertex] = new Path(vertex, -1, Infinity);
				}
			}
			
			//节点是否在队列中
			boolean[] inQueue = new boolean[diNumOfVertice];
			Queue<Integer> queue = new LinkedList<>();
			queue.add(superSource);
			inQueue[superSource] = true;
			
			while (!queue.isEmpty()) {
				int start = queue.poll();
				inQueue[start] = false;
				
				for (Arc arc : graphDi.arcs.get(start)) {
					int end = arc.end;
					//start,end能够到达起点
					if (paths[start].pre != -1 && remain[start][end] > 0) {
						//可以更新end的距离
						if (paths[start].distance + arc.cost < paths[end].distance) {
							paths[end].distance = paths[start].distance + arc.cost;
							paths[end].pre = start;
							if (!inQueue[end]) {
								queue.add(end);
								inQueue[end] = true;
							}
						}
					}
				}
			}
			
			int cur = superSink;
			while (cur >= 0) {
				shortestPath.add(cur);
				cur = paths[cur].pre;
			}
			
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
			
			//输出有向图最短路径
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
		
		//是否满足消费节点的需求
		boolean enough = true;
		for (Consumer consumer : consumers) {
			int sink = consumer.link + numOfVertice;
			if (flowDi[sink][superSink] != consumer.demand) {
				enough = false;
			}
			
				//输出每个消费节点的需求和实际流量
//				writer.println("costumerout "+sink+" need "+consumer.demand+
//						"have " + flowDi[sink][superSink]);
		}
		
		if (enough) {
			int sumCost = sources.size() * serverCost;
			
			for (int i = numOfVertice; i < superSource; i++) {
				for (int j = 0; j < i - numOfVertice; j++) {
					sumCost += (flowDi[i][j] + flowDi[j + numOfVertice][i - numOfVertice]) 
							* cost[i][j];
				}
			}
			
			if (sumCost < minCost) {
				minCost = sumCost;
				minFlow = flowDi;
				minSource = sources;
			}
		}
		
		//删除有向图超级源点和源点的弧及其反向弧
		arcsDi.remove(arcsDi.size() - 2);
		arcsDi.add(arcsDi.size() - 1, new ArrayList<Arc>(defaultSize));
		
		for (int inSource : sources) {
			arcsDi.get(inSource).remove(arcsDi.get(inSource).size() - 1);
		}
	}

	/**
	 * 用最短路径增广最小费用流
	 * @param graphDi
	 * @param remain
	 * @return
	 */
	@SuppressWarnings("unused")
	private static Path[] bellmanFord(GraphDi graphDi, int[][] remain) {
		Path[] paths = new Path[diNumOfVertice];
		for (Arc arc : graphDi.arcs.get(superSource)) {
			paths[arc.end] = new Path(arc.end, superSource, arc.cost);
		}
		for (int vertex = 0; vertex < paths.length; vertex++) {
			if (vertex == superSource) {
				paths[vertex] = new Path(vertex, -2, 0);
			} else if (paths[vertex] == null) {
				paths[vertex] = new Path(vertex, -1, Infinity);
			}
		}
		
		//松弛操作
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
		return paths;
	}
	
	/**
	 * k-means聚类寻找源点组合
	 * @param distance
	 * @param situation
	 */
    
	private static List<List<Integer>> getSourceSituation(int[][] distance) {
		//源点组合
    	List<List<Integer>> situation = new ArrayList<>(numOfConsumers);
    	
    	//开始聚类时存在num个源点
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
			
			//寻找每个汇点最近的源点
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
					map.put(minSource, new ArrayList<Integer>(defaultSize));
				}
				map.get(minSource).add(consumerId);
			}
			
			boolean move = true;
			while (move) {
				move = false;
				Map<Integer, List<Integer>> newMap = new HashMap<>(num);
				
				//对每一组汇点选择新的中心作为新的源点
				for (int source : map.keySet()) {
					List<Integer> consumerCollection = map.get(source);
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
					
					//判断是否与其他组的源点相同
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
    	
    	return situation;
	}

	/**
	 * 求原始网络中汇点到顶点的最短距离
	 * @param edges
	 * @param distance
	 */
	private static void shortestDistance(GraphO graphO, int[][] distance) {
		List<List<GraphO.Edge>> edges = graphO.edges;
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
				int dis = getCost(v0, vertex, edges);
				if (dis == Infinity) {
					paths[vertex] = new Path(vertex, -1, dis);
				} else {
					paths[vertex] = new Path(vertex, v0, dis);
				}
			}
		}
		
		//存储没有确定最短路径的节点
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
				if (cur.distance + edge.cost < paths[next].distance) {
					heap.remove(paths[next]);
					paths[next].pre = start;
					paths[next].distance = cur.distance + edge.cost;
					heap.add(paths[next]);
				}
			}
		}
		
		return paths;
	}

	/**
	 * 获取边的费用
	 * @param start
	 * @param end
	 * @param edges
	 * @return
	 */
	private static int getCost(int start, int end, List<List<GraphO.Edge>> edges) {
		for (GraphO.Edge edge : edges.get(start)) {
			if (edge.end == end) {
				return edge.cost;
			}
		}
		return Infinity;
	}
	
    private static GraphO createGraphO(String[] graphContent) {
		String[] tokens = graphContent[0].split(" ");
    	int numOfVertice = Integer.parseInt(tokens[0]);
    	int numOfEdges = Integer.parseInt(tokens[1]);
    	int numOfConsumers = Integer.parseInt(tokens[2]);
    	int serverCost = Integer.parseInt(graphContent[2]);
    	
    	Deploy.numOfVertice = numOfVertice;
		Deploy.diNumOfVertice = numOfVertice * 2 + 2;
    	Deploy.numOfConsumers = numOfConsumers;
    	Deploy.superSource = 2 * numOfVertice;
		Deploy.superSink = superSource + 1;
		
    	Deploy.serverCost = serverCost;
    	Deploy.numOfConsumers = numOfConsumers;
    	
    	Deploy.defaultSize = (int) (Math.log(numOfConsumers));
    	
    	GraphO graphO = new GraphO(numOfEdges);
    	List<List<GraphO.Edge>> edges = new ArrayList<>(numOfVertice);
    	graphO.edges = edges;
    	for (int i = 0; i < numOfVertice; i++) {
			edges.add(new ArrayList<GraphO.Edge>(defaultSize));
		}
    	
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
    	
		//有向图中除了超级源汇点之外的弧的费用
		cost = new int[2 * numOfVertice][2 * numOfVertice];
		for (List<GraphO.Edge> links : graphO.edges) {
			for (GraphO.Edge edge : links) {
				cost[edge.start + numOfVertice][edge.end] = edge.cost;
			}
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
			if (arcsDi == null) {
				arcsDi = new ArrayList<>(diNumOfVertice);
				
				for (int i = 0; i < diNumOfVertice; i++) {
					arcsDi.add(new ArrayList<Arc>(defaultSize));
				}
				
				//增加入节点到出节点的弧
				for (int vertex = 0; vertex < numOfVertice; vertex++) {
					arcsDi.get(vertex).add(new Arc(vertex, vertex + numOfVertice, Infinity, 0));
				}
				
				//无向边拆为两条有向弧
				for (int start = 0; start < numOfVertice; start++) {
					for (GraphO.Edge edge : graphO.edges.get(start)) {
						int next = edge.end;
						//初始网络中一条边存储为起点和终点相反的两条弧
						if (start < next) {
							int in1 = start;
							int out1 = start + numOfVertice;
							int in2 = next;
							int out2 = next + numOfVertice;
							arcsDi.get(out1).add(new Arc(out1, in2, edge.capacity, edge.cost));
							arcsDi.get(out2).add(new Arc(out2, in1, edge.capacity, edge.cost));
						}
					}
				}
				
				//增加超级汇点和连接汇点的弧
				for (Consumer consumer : consumers) {
					//出节点下标为原始节点下标+原始节点数
					int outSink = consumer.link + numOfVertice;
					arcsDi.get(outSink).add(new Arc(outSink, superSink, consumer.demand, 0));
				}
				
				//添加反向边
				for (List<Arc> links : arcsDi) {
					for (Arc arc : links) {
						if (arc.isback) {
							continue;
						}
						Arc backArc = new Arc(arc.end, arc.start, arc.capacity, -arc.cost);
						backArc.isback = true;
						arcsDi.get(backArc.start).add(backArc);
					}
				}
			}
			
			this.arcs = Deploy.arcsDi;
			
			//增加超级源点到源点的弧,及其反向边
			List<Arc> superSourceLinks = arcs.get(superSource);
			for (int inSource : sources) {
				//入节点下标为原始节点下标
				superSourceLinks.add(new Arc(superSource, inSource, Infinity, 0));
				
				Arc arc = new Arc(inSource, superSource, Infinity, 0);
				arc.isback = true;
				arcs.get(inSource).add(arc);
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
			return "[id=" + id + " link=" + link + " demand=" + demand + "]";
		}
	}

	private static class Path {
		
		private int id;
		
		private int pre;
		
		private int distance;
		
		/**
		 * @param id
		 * @param pre -2:当前节点为最短路径的起点,-1:当前节点未找到连接起点的路径
		 * @param distance
		 */
		public Path(int id, int pre, int distance) {
			this.id = id;
			this.pre = pre;
			this.distance = distance;
		}
		@Override
		public String toString() {
			return "[" + id + "=" + id +" pre= " + pre + " distance= " + distance + "]";
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
