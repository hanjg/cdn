package program325;

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
import java.util.Random;

import com.filetool.util.LogUtil;


public class Deploy
{	
	private static PrintWriter logger;
	
	static {
		try {
			logger = new PrintWriter(new File("data/log.txt"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private static final int Infinity = Integer.MAX_VALUE/2;
	
	/**
	 * (int) (Math.log(numOfConsumers))
	 */
	private static int defaultSize;
	
	/**
	 * 最终输出的路径数量和路径
	 */
	private static List<String> result;
	
	private static List<Integer> minSource;
	
	private static String minCode = "";
	
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
	
	/**
	 * 有向图中超级源点的下标
	 */
	private static int superSource;
	
	/**
	 * 有向图中超级汇点的下标
	 */
	private static int superSink;
	
	/**
	 * 无向图拆点，拆边，添加超级源点，添加超级汇点和对应的弧，得到的有向图中的弧
	 */
	private static List<List<Arc>> arcsDi;
	
	/**
	 * 计时器，从调用deployServer开始计时
	 */
	private static Time timer;
	
	/**
	 * 源点组合的编码和对应的适应度
	 */
	private static HashMap<String, Integer> history;
	
	/**
	 * 代表没有源点的情况，即所有位置字符均为'0'
	 */
	private static StringBuilder defaultCode;
	
	private static String curCode;
	
	private static int curCost;
	
	/**
	 * 超过该费用的个体适应度为0
	 */
	private static int standardCost;

	private static Random random = new Random();
	
	private static GraphO graphO;

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
   		
    	//初始化原始网络和参数
    	init(graphContent);
		
		//输出原始网络的基本信息
    	if (logger != null) {
    		logger.println("numOfVertice:" + numOfVertice +
    				"\nnumOfEdges:" + graphO.numOfEdges + 
    				"\nnumOfConsumers:" + numOfConsumers);
		}
		
		//多次聚类得到结果
//    	sequentialCluster(maxTime);
    	
		//遗传得到结果
		heredity(maxTime);
		
		//计算极端情况的费用：所有消费节点连接服务器
		List<Integer> extreme = new ArrayList<>(numOfConsumers);
		for (Consumer consumer : consumers) {
			extreme.add(consumer.link);
		}
		getMinFlow(extreme);
    	
//    	List<Integer> sources = new ArrayList<>(Arrays.asList(new Integer[]{7, 17, 23, 30, 45, 54, 56, 58, 65, 69, 71, 86, 89, 95, 96, 97, 100, 115, 124, 126, 128, 130, 139, 147, 153}));
//		getMinFlow(sources);
		
		if (logger != null) {
			logger.println("minCost:" + minCost);
			logger.println("minSource:" + minSource);
	    	logger.println("coreTime:" + timer.getTimeDelay() + " ms");
			
			logger.flush();
		}
		
		LogUtil.printLog("---------------------minCost: " + minCost);
		
    	getResult();
    	
//    	check();
    	
    	return result.toArray(new String[result.size()]);
    }

	private static void check() {
		if (logger == null) {
			return;
		}
    	logger.println("check result:");
    	int[][] flow = new int[numOfVertice][numOfVertice];
    	for (int line = 2; line < result.size(); line++) {
			String path = result.get(line);
			String[] tokens = path.split(" ");
			int amount = Integer.parseInt(tokens[tokens.length - 1]);
			for (int i = 0; i + 1 < tokens.length - 2; i++) {
				int from = Integer.parseInt(tokens[i]);
				int to = Integer.parseInt(tokens[i + 1]);
				flow[from][to] += amount;
			}
		}

    	boolean overflow = false;
		for (int from = 0; from < numOfVertice; from++) {
			for (int to = from; to < numOfVertice; to++) {
				if (flow[from][to] == 0 && flow[to][from] == 0) {
					continue;
				}
				int capacity = 0;
				for (GraphO.Edge edge : graphO.edges.get(from)) {
					if (edge.end == to) {
						capacity = edge.capacity;
						break;
					}
				}
				
				if (flow[from][to] > capacity) {
					overflow = true;
				}
//				logger.println(from + "--" + to + " flow: " + flow[from][to] + " " + flow[to][from] + 
//						" capacity is: " + capacity);
//				logger.println(from + "<->" + to + " flow: " + 
//						minFlow[from + numOfVertice][to] + " " + minFlow[to + numOfVertice][from]);
			}
		}
		logger.flush();
		
		boolean enough = true;
		List<Integer> less = new ArrayList<>();
		for (Consumer consumer : consumers) {
			int link = consumer.link;
			int sum = 0;
			for (int from = 0; from < numOfVertice; from++) {
				sum += flow[from][link];
			}
			if (minSource.contains(link)) {
				sum = consumer.demand;
			}
			if (sum < consumer.demand) {
				enough = false;
				less.add(consumer.id);
			}	
		}
		
		logger.println("overflow " + overflow);
		logger.println("enough " + enough);
		if (!enough) {
			logger.println(less);
		}
		logger.println("consumer " + numOfConsumers +" fact " + getConsumerId.size());
		logger.flush();
	}

	private static void heredity(int maxTime) {
		//参数初始化
//		int numOfPopulation = (((numOfVertice / 2) & 1) == 0) ? 
//				numOfVertice /2 : numOfVertice / 2 +1;
//		int numOfPopulation = ((numOfVertice & 1) == 0) ? numOfVertice : numOfVertice + 1;
//		int numOfPopulation = numOfVertice * 2;
		int numOfPopulation = 50;
		
		String[] population = new String[numOfPopulation];
		int[] fitness = new int[numOfPopulation];
		
		final int stopConst = numOfVertice;
//		final int stopConst = numOfVertice * numOfVertice;
//		final int stopConst = 1000;
		final int stopLine = numOfPopulation / 10;
		int stopCount = 0;
//		final int enhanceLine = 0;
		final int enhanceLine = numOfPopulation * 2 / 5;
		final int maxGenerationConst = 100000;
		
		//消费节点(可以看做和汇点重合)和顶点的有效距离(i为消费节点id，j为顶点下标)
    	int[][] distance = new int[numOfConsumers][numOfVertice];
		shortestDistance(distance);
		
		standardCost = serverCost * numOfConsumers * 2;
		
		//初始化种群直到至少有一个个体适应度非0
		boolean inited = false;
		int initCount = 0;
		int index = 0;
		int finished = 0;
		//初始化[finished,len)的种群
		while (!inited) {
			initCount++;
			
			index = finished;
			
//			LogUtil.printLog("now initCount is " + initCount);
			
			//聚类初始化种群
			List<List<Integer>> situation = null;
			situation = kMeans(distance, numOfPopulation / 2);
			for (List<Integer> source : situation) {
				if (index < numOfPopulation / 2) {
					StringBuilder builder = new StringBuilder(defaultCode);
					for (int gene : source) {
						builder.setCharAt(gene, '1');
					}
					population[index++] = builder.toString();
				} else {
					break;
				}
			}
			//防止聚类全部为0,加入极端直连情况
//			StringBuilder builder = new StringBuilder(defaultCode);
//			for (Consumer consumer : consumers) {
//				builder.setCharAt(consumer.link, '1');
//			}
//			population[numOfPopulation - 1] = builder.toString();

			//从极端情况移动源点初始化种群
			List<Integer> sources = new ArrayList<>(defaultSize);
			StringBuilder builder2 = new StringBuilder(defaultCode);
			boolean include = true;
			for (Consumer consumer : consumers) {
				if (include) {
					builder2.setCharAt(consumer.link, '1');
					sources.add(consumer.link);
				} 
				include = !include;
			}
			population[index++] = builder2.toString();
			
			while (index < numOfPopulation) {
				//随机移动源点至相邻节点
				char[] genes = population[index - 1].toCharArray();
				int fromIndex = random.nextInt(sources.size());
				int from = sources.get(fromIndex);
				List<GraphO.Edge> links = graphO.edges.get(from);
				
				int toIndex = 0;
				int to = 0;
				int count = 0;
				while (true) {
					toIndex = random.nextInt(links.size());
					to = links.get(toIndex).end;
					if (!sources.contains(to) || count++ == 5) {
						break;
					}
				}
			
				genes[from] = '0';
				genes[to] = '1';
				population[index++] = new String(genes);
				
				sources.remove(fromIndex);
				sources.add(to);
			}
			
			//随机初始化种群
//			while (index < numOfPopulation) {
//				int numOfSources = random.nextInt(numOfConsumers) + 1;
//				StringBuilder code = new StringBuilder(defaultCode);
//				while (numOfSources > 0) {
//					int pos = random.nextInt(numOfVertice) ;
//					if (code.charAt(pos) == '0') {
//						code.setCharAt(pos, '1');
//						numOfSources--;
//					}
//				}
//				population[index++] = code.toString();
//			}

			//计算各个个体的适应度
			fitness = new int[numOfPopulation];
			for (int i = finished; i < numOfPopulation; i++) {
				
//				LogUtil.printLog("cal fintness " + i + "code is " + population[i]);
				
				curCode = population[i];
				if (history.containsKey(curCode)) {
					fitness[i] = history.get(curCode);
				} else {
					getNewFitness(fitness, i);
				}
			}
			
			LogUtil.printLog("finish cal fitness ");
			
			//判断是否完成初始化
//			for (int i = finished; i < numOfPopulation; i++) {
//				if (fitness[i] != 0) {
//					population[finished++] = population[i];
//				}
//			}
//			if (finished > 0) {
//				break;
//			}
			
			for (int fit : fitness) {
				if (fit != 0) {
					inited = true;
					break;
				}
			}
		}		
		
//		printOriginalGeneration(population);
	
		LogUtil.printLog("----------------------initTimes: " + initCount);

		//遗传迭代
		int generation = 0;
		for (generation = 0; generation < maxGenerationConst; generation++) {
			if (timer.getTimeDelay() > maxTime) {
				break;
			}
			
//			LogUtil.printLog("now is generation " + generation);
			
			//计算各个个体的适应度
			fitness = new int[numOfPopulation];
			//这一代最差的个体
			int minIndex = 0;
			int minFitness = Infinity;
			for (int i = 0; i < numOfPopulation; i++) {
				curCode = population[i];
				if (history.containsKey(curCode)) {
					fitness[i] = history.get(curCode);
				} else {
					getNewFitness(fitness, i);
				}
				
				if (fitness[i] < minFitness) {
					minIndex = i;
					minFitness = fitness[minIndex];
				}
			}
			
			//检查是否收敛
			//最优个体的数量
			int minCount = 0;
			for (String unit : population) {
				if (unit.equals(minCode)) {
					minCount++;
				}
			}
			
			if (minCount > stopLine) {
				stopCount++;
			} else {
				stopCount = 0;
			}
			
			//最优个体是否变化
//			if (minCost == preMinCost) {
//				stopCount++;
//			} else {
//				stopCount = 0;
//			}
			
//			if (stopCount == stopConst / 2) {
//				LogUtil.printLog("---------------------stopCount:" + stopCount + " const:" + stopConst);
//			}
			
			//是否达到终止条件
			if (stopCount > stopConst) {
				break;
			}
			
			
			if (minCount < enhanceLine) {
				//历史最优个体替换当代最差个体
				population[minIndex] = minCode;
				fitness[minIndex] = history.get(minCode);
				//历史最优个体替换当代随机个体
//				population[random.nextInt(numOfPopulation)] = minCode;
			}

			//适应度尺度变换
			//增大适应度最大的个体的影响力
			if (minCount < enhanceLine) {
				int historyMinFitness = history.get(minCode);
				for (int i = 0; i < population.length; i++) {
					if (fitness[i] == historyMinFitness) {
						fitness[i] = fitness[i] * 2;
					}
				}
			}
			//增强适应度大的个体，减弱适应度小的个体
			if (minCount < enhanceLine) {
				int maxFitness = history.get(minCode);
				if (maxFitness > minFitness) {
					double slope = 5.0 * maxFitness / (maxFitness - minFitness);
					for (int i = 0; i < numOfPopulation; i++) {
						fitness[i] = (int) (slope * (fitness[i] - minFitness) + minFitness / 2);
					}
				}
			}
			
			//输出种群信息
//			if (generation < 3000) {
//				printpopulation(population, generation, true);
//				logger.flush();
//			}
			
			
			String[] nextGeneration = new String[numOfPopulation];
			
			//选择
			//[0,i)的适应度之和
			int[] sumFitness = new int[numOfPopulation + 1];
			for (int i = 1; i <sumFitness.length; i++) {
				sumFitness[i] += sumFitness[i - 1] + fitness[i - 1];
			}
			int totalFitness = sumFitness[numOfPopulation];
			
			//判断totalfitness是否溢出
			if (totalFitness < 0) {
				LogUtil.printLog("totalFitness :" + totalFitness);
			}
			
			//赌盘选择numOfPopulation个体
			for (int i = 0; i < numOfPopulation; i++) {
				//选择i则sumfitness[i]<=ran<sumfitness[i+1]
				int ran = random.nextInt(totalFitness);
				int sel = 0;
				//ran<sumFitness[numOfPopulation]=totalFitness
				while (ran >= sumFitness[sel + 1]) {
					sel++;
				}
				nextGeneration[i] = population[sel];
			}
			
//			double crossProbability = minCount < enhanceLine ? 0.75 : 1;
			double crossProbability = (0.5 / numOfPopulation) * minCount + 0.5;
			
			//分组交叉
			for (int i = 0; i < numOfPopulation; i += 2) {
				if (Math.random() > crossProbability) {
					continue;
				}
				
//				if (nextGeneration[i].equals(nextGeneration[i + 1])) {
//					continue;
//				}
				
				//单点交叉cut,cut+1之间cut[0,len-1);
				int cut = random.nextInt(nextGeneration[0].length() + 1);
				String temp = nextGeneration[i].substring(cut);
				nextGeneration[i] = nextGeneration[i].substring(0, cut) + 
						nextGeneration[i + 1].substring(cut);
				nextGeneration[i + 1] = nextGeneration[i + 1].substring(0, cut) + temp;
				
				//之前错误，交换时未用temp。。
//				nextGeneration[i] = nextGeneration[i].substring(0, cut1 + 1) 
//						+ nextGeneration[i + 1].substring(cut1 + 1, nextGeneration[i].length());
//				nextGeneration[i + 1] = nextGeneration[i + 1].substring(0, cut1 + 1)
//						+ nextGeneration[i].substring(cut1 + 1, nextGeneration[i].length());
				
				//两点交叉
//				int temp1 = random.nextInt(nextGeneration[0].length() + 1);
//				int temp2 = random.nextInt(nextGeneration[0].length() + 1);
//				int cut1 = Math.min(temp1, temp2);
//				int cut2 = Math.max(temp1, temp2);
//
//				String mid = nextGeneration[i].substring(cut1, cut2);
//				
//				nextGeneration[i] = nextGeneration[i].substring(0, cut1) +
//						nextGeneration[i + 1].substring(cut1, cut2) + 
//						nextGeneration[i].substring(cut2);
//				
//				nextGeneration[i + 1] = nextGeneration[i + 1].substring(0, cut1) +
//						mid + nextGeneration[i + 1].substring(cut2);
				
				//均匀交叉
//				StringBuilder builder1 = new StringBuilder(nextGeneration[0].length());
//				StringBuilder builder2 = new StringBuilder(nextGeneration[0].length());
//				char[] ca1 = nextGeneration[i].toCharArray();
//				char[] ca2 = nextGeneration[i + 1].toCharArray();
//				for (int j = 0; j < ca1.length; j++) {
//					if (Math.random() < 0.5) {
//						builder1.append(ca1[j]);
//						builder2.append(ca2[j]);
//					} else {
//						builder1.append(ca2[j]);
//						builder2.append(ca1[j]);
//					}
//				}
//				
//				nextGeneration[i] = builder1.toString();
//				nextGeneration[i + 1] = builder2.toString();
			}
			
			//单亲交叉
//			for (int i = 0; i < numOfPopulation; i++) {
//				if (Math.random() > crossProbability) {
//					continue;
//				}
//				
//				int temp1 = random.nextInt(nextGeneration[0].length() + 1);
//				int temp2 = random.nextInt(nextGeneration[0].length() + 1);
//				int cut1 = Math.min(temp1, temp2);
//				int cut2 = Math.max(temp1, temp2);
//				int len = random.nextInt(Math.min(cut2 - cut1, numOfPopulation - cut2) + 1);
//				String string = population[i];
//				population[i] = string.substring(0, cut1) + string.substring(cut2, cut2 + len) +
//						string.substring(cut1 + len, cut2) + string.substring(cut1, cut1 + len) +
//						string.substring(cut2 + len);
//			}
			
			//变异
//			double mutateProbability = minCount < enhanceLine ? 0.05 : 0.2;
			double mutateProbability = (0.2 / numOfPopulation) * minCount + 0.05;
			
			for (int i= 0; i < numOfPopulation; i++) {
				char[] genes = nextGeneration[i].toCharArray();
				if (Math.random() < mutateProbability) {
					if (population[i].equals(minCode)) {
						//随机移动源点至相邻节点
						int from = minSource.get(random.nextInt(minSource.size()));
						List<GraphO.Edge> links = graphO.edges.get(from);
						int to = links.get(random.nextInt(links.size())).end;
						
						genes[from] = '0';
						genes[to] = '1';
						nextGeneration[i] = new String(genes);
					} else {
						//随机增删源点
						int pos = random.nextInt(nextGeneration[0].length());
						genes[pos] = genes[pos] == '0' ? '1' : '0';
						nextGeneration[i] = new String(genes);
					}
					
					//随机增删源点
//					int pos = random.nextInt(nextGeneration[0].length());
//					genes[pos] = genes[pos] == '0' ? '1' : '0';
//					nextGeneration[i] = new String(genes);
				
				} 
			}
			
			population = nextGeneration;
			
		
		}
		
		if (logger != null) {
			logger.println("generationCount :" + generation);
		}
		LogUtil.printLog("-----------------------generationCount: " + generation);
		
	}

	/**
	 * @param population
	 */
	private static void printOriginalGeneration(String[] population) {
		if (logger == null) {
			return;
		}
		List<Integer> sources;
		logger.println("------------------------");
		logger.println("original generation:");
		for (String unit : population) {
			sources = new ArrayList<>(defaultSize);
			char[] genes = unit.toCharArray();
			for (int gene = 0; gene < genes.length; gene++) {
				if (genes[gene] == '1') {
					sources.add(gene);
				}
			}
			logger.println("fitness : " + history.get(unit) + "  " + sources);
		}
		logger.println("-------------------------------");
		logger.flush();
	}

	/**
	 * 获取之前未出现过的个体的适应度
	 * @param fitness
	 * @param position
	 */
	private static void getNewFitness(int[] fitness, int position) {
		List<Integer> sources = new ArrayList<>(defaultSize);
		char[] genes = curCode.toCharArray();
		for (int gene = 0; gene < genes.length; gene++) {
			if (genes[gene] == '1') {
				sources.add(gene);
			}
		}
		
		getMinFlow(sources);
		
		int fit = curCost > standardCost ? 0 : (standardCost - curCost) / 2;
		
		fitness[position] = fit;
		history.put(curCode, fit);
	}

	/**
	 * @param population
	 * @param count
	 */
	private static void printpopulation(String[] population, int count, boolean detail) {
		if (logger == null) {
			return;
		}
		logger.println("generation " + count + " :");

		int minCount = 0;
		for (String unit : population) {
			if (unit.equals(minCode)) {
				minCount++;
			}
		}
		logger.println("minCost :" + minCost + " num :" + minCount
				+ " maxFitness :" + history.get(minCode) );
		
		if (detail) {
//			for (String unit : population) {
//				List<Integer> sources = new ArrayList<>(defaultSize);
//				char[] genes = unit.toCharArray();
//				for (int gene = 0; gene < genes.length; gene++) {
//					if (genes[gene] == '1') {
//						sources.add(gene);
//					}
//				}
//				logger.println("fitness : " + history.get(unit) + "  " + sources);
//			}
			List<Integer> amount = new ArrayList<>(defaultSize);
			//适应度和amount中的下标
			Map<Integer, Integer> map = new HashMap<>(defaultSize);
			for (String unit : population) {
				int fit = history.get(unit);
				if (map.containsKey(fit)) {
					amount.set(map.get(fit), amount.get(map.get(fit)) + 1);
				} else {
					map.put(fit, amount.size());
					amount.add(1);
				}
			}
			logger.println("varity : " + map.size());
			for (int key : map.keySet()) {
				logger.print(key + ":" + amount.get(map.get(key)) + " ");
			}
			logger.println();
		}
		
		logger.println("---------------------------------");
	}

	private static void sequentialCluster(int maxTime) {
		//消费节点(可以看做和汇点重合)和顶点的有效距离(i为消费节点id，j为顶点下标)
    	int[][] distance = new int[numOfVertice][numOfVertice];
		shortestDistance(distance);
		
		final int maxClusterConst = 100000;
//		final int stopConst = 10 * numOfVertice * numOfVertice;
		final int stopConst = maxClusterConst;
		int stopCount = 0;
		
		int count = 0;
		for (count = 0; count < maxClusterConst; count++) {
			if (timer.getTimeDelay() > maxTime) {
				break;
			}
			
			//前一代最优解
			int originMinCost = minCost;
			
			//获取源点的组合
			List<List<Integer>> situation = kMeans(distance, numOfConsumers);
			
			//计算该源点组合下的最小费用流
	    	for (List<Integer> sources : situation) {
	    		
	    		//源点组合的编码
				StringBuilder builder = new StringBuilder(defaultCode);
				for (int source : sources) {
					builder.setCharAt(source, '1');
				}
				curCode = builder.toString();
				
				//编码出现过则跳过计算最小费用流的步骤
				if (history.containsKey(curCode)) {
					continue;
				}
				
	    		getMinFlow(sources);
	    		
	    		//记录计算过最小费用流的编码
	    		history.put(curCode, 0);
	    	}
	    	
	    	//判断是否收敛
	    	if (minCost == originMinCost) {
				stopCount++;
			} else {
				stopCount = 0;
			}
	    	
	    	if (stopCount > stopConst) {
				break;
			}
	    }
		
		LogUtil.printLog("\nclusterCount: " + count);
	}
    
	/**
	 * 通过最小费用流矩阵求得路径
	 * @param graphO
	 * @param edges
	 */
	private static void getResult() {
		List<List<GraphO.Edge>> edges = graphO.edges;
		
		//修改原始网络，增加虚拟源汇点和对应的边(在原始网络上修改)
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
		
    	//原始网络中添加超级源汇点，得到新的流网络
    	//复制原始网络的边
    	int[][] remain = new int[numOfVertice + 2][numOfVertice + 2];
    	for (int i = 0; i < numOfVertice; i++) {
			for (int j =0; j < numOfVertice; j++) {
				remain[i][j] = flow[i][j];
			}
		}
    	
    	//在流网络中增加超级源汇点对应的边,保证源点的流量足够,汇点流量为需求
    	for (int source : minSource) {
			remain[numOfVertice][source] = Infinity;
		}
    	for (Consumer consumer : consumers) {
			remain[consumer.link][numOfVertice + 1] = consumer.demand;
		}
    	
    	result = new ArrayList<>(defaultSize);

    	while (true) {
			if (!dfs(new boolean[numOfVertice + 2], numOfVertice, Infinity, new ArrayList<Integer>(defaultSize), remain)) {
				break;
			}
		}
    	
    	int count = result.size();
    	
    	result.addAll(0, Arrays.asList(new String[]{Integer.toString(count),""}));
	}
	
	private static void printMinFlow() {
		if (logger == null) {
			return;
		}
		
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
     * @param visited
     * @param cur
     * @param flow 源点到cur的路径流量
     * @param path
     * @return 是否找到虚拟源点到虚拟汇点的路径
     */
    private static boolean dfs(boolean[] visited, int cur, int flow, List<Integer> path, int[][] remain) {
    	if (cur == numOfVertice + 1) {
    		path.add(cur);
			StringBuilder builder = new StringBuilder();
			for (int i = 1; i < path.size() - 1; i++) {
				builder.append(path.get(i) + " ");
			}
			int consumerId = getConsumerId.get(path.get(path.size() - 2));
			builder.append(consumerId + " ");
			builder.append(flow);
			
			result.add(builder.toString());
			
			for (int i = 0; i < path.size() - 1; i++) {
				remain[path.get(i)][path.get(i+1)] -= flow;
//				expandFlow[path.get(i)][path.get(i+1)] += flow;
//				int from = path.get(i);
//				int to = path.get(i+1);
//				int capacity = 0;
//				for (GraphO.Edge edge : graphO.edges.get(from)) {
//					if (edge.end == to) {
//						capacity = edge.capacity;
//						break;
//					}
//				}
//				logger.println("expandFlow " + from + "->" + to + " :" + expandFlow[from][to] + " capacity: " + capacity);
			}
			
			
			return true;
		}
    	
    	path.add(cur);
    	visited[cur] = true;
    	
    	for (GraphO.Edge edge : graphO.edges.get(cur)) {
			int next = edge.end;
			if (!visited[next] && remain[cur][next] > 0) {
				if (dfs(visited, next, Math.min(flow, remain[cur][next]), path, remain)) {
					return true;
				}
			}
		}
    	path.remove(path.size() - 1);
    	visited[cur] = false;
    	
    	return false;
    }
    
    private static boolean dfs2(boolean[] visited, int cur, int flow, List<Integer> path, int[][] remain) {
    	if (cur == numOfVertice + 1) {
    		path.add(cur);
			StringBuilder builder = new StringBuilder();
			for (int i = 1; i < path.size() - 1; i++) {
				builder.append(path.get(i) + " ");
			}
			int consumerId = getConsumerId.get(path.get(path.size() - 2));
			builder.append(consumerId + " ");
			builder.append(flow);
			
			result.add(builder.toString());
			
			for (int i = 0; i < path.size() - 1; i++) {
				remain[path.get(i)][path.get(i+1)] -= flow;
			}
			return true;
		}
    	
    	path.add(cur);
    	visited[cur] = true;
    	
    	int next = Infinity;
    	int max = Infinity;
    	while (true) {
    		next = -1;
    		max = -1;
    		for (GraphO.Edge edge : graphO.edges.get(cur)) {
				int end = edge.end;
				if (!visited[end] && remain[cur][end] > 0) {
					if (remain[cur][end] > max) {
						next = end;
						max = remain[cur][end];
					}
				}
			}
    		if (next == -1) {
				break;
			} else {
				if (dfs(visited, next, Math.min(flow, remain[cur][next]), path, remain)) {
					return true;
				}
			}
		}
		
    	path.remove(path.size() - 1);
    	visited[cur] = false;
    	
    	return false;
    } 
    
	/**
	 * 根据源点组合建立有向图，求最小费用流
	 * @param graphO
	 * @param result
	 * @param sources
	 */
	private static void getMinFlow(List<Integer> sources) {
		//根据初始无向网络和源点位置，拆点，拆边，增加超级源汇点和对应弧，新建有向图
		GraphDi graphDi = new GraphDi(sources);
		
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
			
			//无最短路径，已获得最小费用流
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
//			if (logger != null) {
//				logger.println("flow:"+minFlow);
//				logger.print("shortestPath:");
//				for (int i = shortestPath.size() - 1; i >= 0; i--) {
//					logger.print(shortestPath.get(i)+" ");
//				}
//				logger.println();
//			}
			
			//求剩余网络
			for (int i = shortestPath.size() - 1; i >= 1; i--) {
				int start = shortestPath.get(i);
				int end = shortestPath.get(i-1);
				Arc section = null;
				for (Arc arc : graphDi.arcs.get(start)) {
//					if (arc.start == start) {错误
					if (arc.end == end) {
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
					//反向弧start->end减少正向弧end-start的流量(bug已修正)
					flowDi[end][start] -= minFlow;
					remain[end][start] += minFlow;
					remain[start][end] = flowDi[end][start];
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

//			if (consumer.id == 0) {
//				logger.println(consumer);
//				logger.println("flowDi[" + sink + "][" + superSink + "]=" + flowDi[sink][superSink]);
//			}
			//输出每个消费节点的需求和实际流量
//			if (logger != null) {
//				logger.println("costumerout "+sink+" need "+consumer.demand+
//						"have " + flowDi[sink][superSink]);
//			}
		}
		
		int sumCost = Infinity;
		if (enough) {
			sumCost = sources.size() * serverCost;
			
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
				minCode = curCode;
			}
		}

		curCost = sumCost;
		
		//删除有向图超级源点到源点的弧及其反向弧
		arcsDi.remove(arcsDi.size() - 2);
		arcsDi.add(arcsDi.size() - 1, new ArrayList<Arc>(defaultSize));
		
		for (int inSource : sources) {
			arcsDi.get(inSource).remove(arcsDi.get(inSource).size() - 1);
		}
		
	}

	/**
	 * bellmanFord最短路径增广最小费用流
	 * @param graphDi
	 * @param remain
	 * @return
	 */
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
	 * @param need 需要的源点组合个数
	 * @return
	 */
	private static List<List<Integer>> kMeans(int[][] distance, int need) {
		//源点组合
    	List<List<Integer>> situation = new ArrayList<>(numOfConsumers);
    	
    	//开始聚类时存在num个源点
    	while (need-- > 0) {
			int num = random.nextInt(numOfConsumers - 1) + 1;

    		//源点的位置
			List<Integer> newSources = new ArrayList<>(num);
			while (newSources.size() < num) {
				int temp = random.nextInt(numOfVertice);
				if (!newSources.contains(temp)) {
					newSources.add(temp);
				}
			}
			//K=源点下标 V=属于该源点的消费节点id的集合
			Map<Integer, List<Integer>> map = new HashMap<>(num);
			
			//寻找距离汇点最近的源点，则该汇点属于该源点
			for (int consumerId = 0; consumerId < numOfConsumers; consumerId++) {
				int minSource = -1;
				int minDistance = Infinity;
				for (int source : newSources) {
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
			
			//标准k-means，效果不好，原因未知。。。
//			while (move) {
//
////				if (logger != null) {
////					logger.println(map.keySet());
////				}
//				
//				move = false;
//				
//				//每一组汇点寻找新的中心作为新的源点
//				newSources = new ArrayList<>(num);
//				for (int source : map.keySet()) {
//					int newCenter = -1;
//					int minSum = Infinity;
//					for (int vertex = 0; vertex < numOfVertice; vertex++) {
//						int sum = 0;
//						for (int consumerId : map.get(source)) {
//							sum += distance[consumerId][vertex];
//						}
//						if (sum < minSum) {
//							newCenter = vertex;
//							minSum = sum;
//						}
//					}
//					
//					if (!newSources.contains(newCenter)) {
//						newSources.add(newCenter);
//					}
//					
//					if (source != newCenter) {
//						move = true;
//					} 
//				}
//				
//				map = new HashMap<>(num);
//				
//				//寻找距离汇点最近的源点，则该汇点属于该源点
//				for (int consumerId = 0; consumerId < numOfConsumers; consumerId++) {
//					int minSource = -1;
//					int minDistance = Infinity;
//					for (int source : newSources) {
//						if (distance[consumerId][source] < minDistance) {
//							minSource = source;
//							minDistance = distance[consumerId][minSource];
//						}
//					}
//					if (!map.containsKey(minSource)) {
//						map.put(minSource, new ArrayList<Integer>(defaultSize));
//					}
//					map.get(minSource).add(consumerId);
//				}
//			}
			
//			if (logger != null) {
//				logger.println(map.keySet());
//				logger.println("----------------------------------------");
//			}
			
			
			//改编的k-means
			//初始化在一个集合的汇点直到最终也在同一个集合
			while (move) {
				move = false;
				Map<Integer, List<Integer>> newMap = new HashMap<>(num);
				
				//对每一组汇点选择新的中心作为新的源点
				for (int source : map.keySet()) {
					List<Integer> sinkCollection = map.get(source);
					int newCenter = -1;
					int minSum = Infinity;
					for (int vertex = 0; vertex < numOfVertice; vertex++) {
						int sum = 0;
						for (int consumerId : sinkCollection) {
							sum += distance[consumerId][vertex];
						}
						if (sum < minSum) {
							newCenter = vertex;
							minSum = sum;
						}
					}
					if (source != newCenter) {
						move = true;
					} 
					
					//判断是否与其他组的源点相同
					if (newMap.containsKey(newCenter)) {
						newMap.get(newCenter).addAll(sinkCollection);
					} else {
						newMap.put(newCenter, sinkCollection);
					}
				}
				map = newMap;
			}
			
			situation.add(new ArrayList<>(map.keySet()));
		
		}
    	
    	return situation;
	}

	/**
	 * 求原始网络中汇点(和消费节点相连的节点)到顶点的最短距离
	 * @param edges
	 * @param distance
	 */
	private static void shortestDistance(int[][] distance) {
		List<List<GraphO.Edge>> edges = graphO.edges;
		for (int id = 0; id < numOfConsumers; id++) {
    		Consumer consumer = consumers.get(id);
			int sink = consumer.link;
			
			Path[] paths = dijstra(edges, sink);
			
			for (int i = 0; i < numOfVertice; i++) {
//				distance[id][i] = (paths[i].distance << 16) / consumer.demand;
				distance[id][i] = paths[i].distance;
			}
		}
	}

	/**
	 * heap + dijstra
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
//				int dis = getCost(v0, vertex, edges);
//				if (dis == Infinity) {
//					paths[vertex] = new Path(vertex, -1, dis);
//				} else {
//					paths[vertex] = new Path(vertex, v0, dis);
//				}
				paths[vertex] = new Path(vertex, -1, getEdgeCost(v0, vertex, edges));
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
			Path path = heap.poll();
			int start = path.id;
			for (GraphO.Edge edge : edges.get(start)) {
				int end = edge.end;
				if (path.distance + edge.cost < paths[end].distance) {
					heap.remove(paths[end]);
					paths[end].pre = start;
					paths[end].distance = path.distance + edge.cost;
					heap.add(paths[end]);
				}
			}
		}
		
		return paths;
	}

	private static int getEdgeCost(int start, int end, List<List<GraphO.Edge>> edges) {
		for (GraphO.Edge edge : edges.get(start)) {
			if (edge.end == end) {
				return edge.cost;
			}
		}
		return Infinity;
	}
	
    /**
     * 构造原始无向网络，并初始化参数
     * @param graphContent
     * @return
     */
    private static void init(String[] graphContent) {
		String[] tokens = graphContent[0].split(" ");
    	int numOfVertice = Integer.parseInt(tokens[0]);
    	int numOfEdges = Integer.parseInt(tokens[1]);
    	int numOfConsumers = Integer.parseInt(tokens[2]);
    	int serverCost = Integer.parseInt(graphContent[2]);
    	
    	Deploy.serverCost = serverCost;
    	Deploy.numOfConsumers = numOfConsumers;
    	Deploy.defaultSize = (int) (Math.log(numOfConsumers));
    	
    	//初始化原始无向图网络
    	Deploy.graphO = new GraphO(numOfEdges);
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
//			for (String string : edge) {
//				System.out.print(string + " ");
//			}
//			System.out.println();
			int v1 = Integer.parseInt(edge[0]);
			int v2 = Integer.parseInt(edge[1]);
			int capacity = Integer.parseInt(edge[2]);
			int cost = Integer.parseInt(edge[3]);
//			//判断是否为重边
//			if (isMutiEdge(v1, v2)) {
//				//通过在重复的边中间添加节点，解决重边问题
//				edges.add(new ArrayList<GraphO.Edge>(defaultSize));
//				int plus = numOfVertice + numOfMultiEdge++;
//				edges.get(v1).add(new GraphO.Edge(v1, plus, capacity, cost));
//				edges.get(plus).add(new GraphO.Edge(plus, v1, capacity, cost));
//				edges.get(plus).add(new GraphO.Edge(plus, v2, Infinity, 0));
//				edges.get(v2).add(new GraphO.Edge(v2, plus, Infinity, 0));
//			}
//			if (isMutiEdge(v1, v2)) {
//				continue;
//			}
			edges.get(v1).add(new GraphO.Edge(v1, v2, capacity, cost));
			edges.get(v2).add(new GraphO.Edge(v2, v1, capacity, cost));
		
		}

    	Deploy.numOfVertice = numOfVertice;
		Deploy.diNumOfVertice = numOfVertice * 2 + 2;
    	Deploy.numOfConsumers = numOfConsumers;
    	Deploy.superSource = 2 * numOfVertice;
		Deploy.superSink = superSource + 1;
    	
    	//初始化消费节点和相关参数
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
    	
		//计算有向图中除了连接超级源汇点之外的弧的费用
		Deploy.cost = new int[2 * numOfVertice][2 * numOfVertice];
		
		for (List<GraphO.Edge> links : graphO.edges) {
			for (GraphO.Edge edge : links) {
				cost[edge.start + numOfVertice][edge.end] = edge.cost;
			}
		}
		
		//初始化map，记录源点组合编码和对应的最小费用
		//初始化默认源点组合的编码
		history = new HashMap<>(numOfVertice * numOfVertice);
		defaultCode = new StringBuilder(numOfVertice);
		for (int i = 0; i < numOfVertice; i++) {
			defaultCode.append(0);
		}
		history.put(defaultCode.toString(), Infinity);
		
	}
    
    private static boolean isMutiEdge(int v1, int v2) {
    	for (GraphO.Edge edge : graphO.edges.get(v1)) {
			if (v2 == edge.end) {
				return true;
			}
		}
    	return false;
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
			
			@Override
			public String toString() {
				return "[" + start + "--" + end + " capacity " + capacity + " cost " + cost + "]";
			}
		}
	}
	private static class GraphDi {
		
		private List<List<Arc>> arcs;
		
		public GraphDi(List<Integer> sources) {
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
