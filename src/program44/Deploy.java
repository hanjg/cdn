package program44;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
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
	 * 最终输出的路径数量和路径
	 */
	private static List<String> result;
	
	private static int[] minSource;
	
	private static int numOfMinSource;
	
	private static String minCode = "";
	
	private static int minCost = Infinity;
	
	/**
	 * 对应边上的流量
	 */
	private static int[] minFlow;
	
    private static int numOfConsumers;

	private static Consumer[] consumers;
	
	/**
	 * 由consumer.link查找consumer.id
	 */
	private static Map<Integer, Integer> getConsumerId;
	
	private static int totalDemand;

	private static int serverCost;

	private static int numOfVertice;
	
	private static int numOfEdges;
	
	private static Arc[] edges;
	
	private static int[] edgeHeadIndex;
	
	private static int[] numOfLinkedEdges;
	
	/**
	 * 包括超级源汇点对应的边和反向边
	 */
	private static int numOfArcs;
	
	/**
	 * 拆点，增加超级源汇点之后的顶点数
	 */
	private static int diNumOfVertice;
	
	private static int superSource;
	
	private static int superSink;
	
	/**
	 * 偶数正向弧，奇数反向弧
	 */
	private static Arc[] arcs;
	
	/**
	 * 无向图中每个顶点的第一个邻接弧的在arcs中的下标<p>
	 * 节点的下标顺序：入节点，出节点，超级源点，超级汇点
	 */
	private static int[] arcHeadIndex;
	
	/**
	 * 源点组合的编码和对应的适应度
	 */
	private static HashMap<String, SourceInfo> history;
	
	/**
	 * 代表没有源点的情况，即所有位置字符均为'0'
	 */
	private static StringBuilder defaultCode;
	
	private static String curCode;
	
	private static int curCost;
	
	private static int[] curSource;
	
	private static int numOfCurSource;
	
	private static int[] curFlow;
	
	/**
	 * 超过该费用的个体适应度为0
	 */
	private static int standardCost;

	private static Random random = new Random();
	
	/**
	 * 计时器，从调用deployServer开始计时
	 */
	private static Time timer = new Time();
	
	private static int maxTime = 80000;
	
	private static int getMinTime = 0;
    
    private static int getMinCount = 0;

	static int Hmove = 0;
	
	static int Hdel = 0;
	
	static int Hadd = 0;
	
	/**
     * 你需要完成的入口
     * <功能详细描述>
     * @param graphContent 用例信息文件
     * @return [参数说明] 输出结果信息
     * @see [类、类#方法、类#成员]
     */
    public static String[] deployServer(String[] graphContent) {
    	//初始化原始网络和参数
    	init(graphContent);
    	
//    	LogUtil.printLog("finish original graph");
    	
		//多次聚类得到结果
//    	sequentialCluster();
    	
		
		//遗传得到结果
		heredity();
		
    	//退火得到结果
//    	coolDown(maxTime);
    	
//    	testGetMinFlow();
		
		//计算极端情况的费用：所有消费节点连接服务器
//		calStraghtLink();
    	
//    	verticeInfo();
    	
    	
		if (logger != null) {
    		logger.println("numOfVertice:" + numOfVertice +
    				"\nnumOfEdges:" + numOfEdges + 
    				"\nnumOfConsumers:" + numOfConsumers);
		
			logger.println("minCost:" + minCost);
			logger.println("numOfMinSources:" + numOfMinSource);
			logger.print("minSource:");
			for (int i = 0; i < numOfMinSource; i++) {
				logger.print(minSource[i] + ",");
			}
			logger.println();
	    	logger.println("coreTime:" + timer.getTimeDelay() + " ms");
	    	if (getMinCount > 0) {
				logger.println("sumTime : " + getMinTime + " count : " + getMinCount + " mean : " + getMinTime / getMinCount);
			}
			logger.flush();
		}
		
		LogUtil.printLog("---------------------minCost: " + minCost +
				"\nsumTime : " + getMinTime + " count : " + getMinCount + " mean : " + getMinTime / getMinCount);

    	result = new ArrayList<>(numOfVertice);
    	getResult();
    	
//    	check();
    	
    	return result.toArray(new String[result.size()]);
    }

	private static List<List<Integer>> levelCluster(int[][] validDistance, int need, int lessThan) {
		
		//汇点之间的有效距离，i->j
		double[][] sinkDistance = new double[numOfConsumers][numOfConsumers];
		for (int id= 0; id < numOfConsumers; id++) {
			for (Consumer consumer : consumers) {
				sinkDistance[consumer.id][id] = validDistance[id][consumer.link];
			}
		}
		
//		int[] centerSinks = new int[numOfConsumers];
//		for (int i = 0; i < numOfConsumers; i++) {
//			centerSinks[i] = -1;
//		}
		//中心汇点-对应汇点组合
		Map<Integer, List<Integer>> map = new HashMap<>(numOfConsumers);
		for (int i = 0; i < numOfConsumers; i++) {
			List<Integer> temp = new ArrayList<>(numOfConsumers);
			temp.add(i);
			map.put(i, temp);
		}
		
    	List<List<Integer>> situation = new ArrayList<>(need);
    	
		while (need > 0) {
			List<Integer> sources = new ArrayList<>(numOfConsumers);
			for (int centerSink : map.keySet()) {
				List<Integer> sinks = map.get(centerSink);
				if (sinks.size() == 1) {
					sources.add(consumers[sinks.get(0)].link);
				} else {
					int minCenter = -1;
					int minSum = Infinity;
					for (int vertex = 0; vertex < numOfVertice; vertex++) {
						int sum = 0;
						for (int sink : sinks) {
							sum += validDistance[sink][vertex];
						}
						if (sum < minSum) {
							minCenter = vertex;
							minSum = sum;
						}
					}
					sources.add(minCenter);
				}
			}
			
			if (sources.size() < lessThan) {
				situation.add(sources);
				need--;
			}
			
			//找出距离最近的两个汇点，合并
			int from = 0;
			int to = 0;
			double min = Double.MAX_VALUE;
			for (int i = 0; i < numOfConsumers; i++) {
				for (int j = 0; j < numOfConsumers; j++) {
					if (sinkDistance[i][j] > 0 && sinkDistance[i][j] < min) {
						from = i;
						to = j;
						min = Math.min(min, sinkDistance[i][j]);
					}
				}
			}
			int comb1 = Math.min(from, to);
			int comb2 = Math.max(from, to);
			
			//comb2与comb1合并
//			centerSinks[comb2] = comb1;
			map.get(comb1).addAll(map.get(comb2));
			map.remove(comb2);
			
			//修正comb1->j的距离
			for (int j = 0; j < numOfConsumers; j++) {
				sinkDistance[comb1][j] = (sinkDistance[comb1][j] + sinkDistance[comb2][j]) / 2;
			}
			//修正i->comb1的距离
			for (int i = 0; i < numOfConsumers; i++) {	
				sinkDistance[i][comb1] = (sinkDistance[i][comb1] + sinkDistance[i][comb2]) / 2;
			}
			//comb2->j的距离参考comb1->j的距离
			for (int j = 0; j < numOfConsumers; j++) {
				sinkDistance[comb2][j] = -1;
			}
			//i->comb2的距离参考i->comb1
			for (int i = 0; i < numOfConsumers; i++) {
				sinkDistance[i][comb2] = -1;
			}
			//comb1和comb2合并，所以距离为0
			sinkDistance[comb1][comb1] = sinkDistance[comb1][comb2] 
					= sinkDistance[comb2][comb1] = sinkDistance[comb2][comb2] = 0;
		}
		
//		for (int i = 0 ; i < situation.size(); i++) {
//			Set<Integer> set = new HashSet<>(situation.get(i));
//			logger.println("listsize " + situation.get(i).size() + situation.get(i));
//			logger.println("setsize " + set.size() + set);
//		}
		
//		for (List<Integer> sources : situation) {
//			curSource = new int[numOfVertice];
//			numOfCurSource = 0;
//			for (int source : sources) {
//				curSource[numOfCurSource++] = source;
//			}
//			
//			StringBuilder builder = new StringBuilder(defaultCode);
//			for (int i = 0; i < numOfCurSource; i++) {
//				builder.setCharAt(curSource[i], '1');
//			}
//			curCode = builder.toString();
//			
//			getMinFlow();
//			
//			LogUtil.printLog("numofsource " + numOfCurSource + " cost " + curCost);
//
//			logger.println("numofsource " + numOfCurSource + " cost " + curCost + sources);
//			logger.flush();
//		}
		
		return situation;
	}

	/**
	 * 计算直连的极端情况
	 */
	private static void calStraghtLink() {
		List<Integer> extreme = new ArrayList<>(numOfConsumers);
		for (Consumer consumer : consumers) {
			extreme.add(consumer.link);
		}
		getMinFlow();
	}

	private static void testGetMinFlow() {
		List<Integer> temp = new ArrayList<>(Arrays.asList(new Integer[]{7, 14, 17, 19, 25, 26, 29, 32, 35, 40, 43, 44, 59, 70, 73, 82, 84, 89, 93, 101, 108, 111, 120, 124, 126, 127, 129, 131, 133, 137, 138, 147, 152, 155, 161, 164, 166, 167, 168, 175, 178, 184, 186, 187, 190, 194, 195, 198, 200, 205, 207, 218, 223, 227, 230, 234, 237, 238, 242, 252, 254, 259, 263, 267, 268, 270, 271, 272, 278, 281, 287, 288, 308, 320, 321, 326, 327, 328, 330, 333, 334, 335, 336, 338, 346, 349, 363, 364, 370, 375, 381, 383, 385, 387, 397, 400, 402, 409, 413, 415, 416, 423, 426, 433, 438, 457, 459, 463, 466, 474, 475, 482, 491, 496, 500, 503, 505, 507, 520, 525, 528, 529, 531, 537, 541, 542, 548, 552, 557, 561, 565, 569, 575, 576, 584, 587, 592, 594, 599, 609, 625, 632, 634, 636, 638, 641, 642, 643, 646, 651, 653, 660, 661, 663, 668, 670, 671, 675, 676, 679, 686, 687, 695, 711, 714, 718, 719, 724, 728, 731, 739, 741, 742, 745, 751, 753, 755, 756, 757, 760, 765, 767, 770, 774, 778, 784, 791, 795, 797}));
//    	List<Integer> temp = new ArrayList<>(Arrays.asList(new Integer[]{166, 320, 62, 463, 604, 413, 551, 53, 349, 520, 438, 684, 16, 707, 423, 200, 241, 290, 146, 82, 291, 646, 739, 383, 548, 89, 145, 346, 331, 641, 791, 760, 742, 769, 594, 731, 482, 500, 267, 503, 719, 19, 440, 43, 625, 660, 270, 576, 299, 178, 408, 653, 370, 584, 670, 402, 259, 753, 529, 765, 334}));
		System.out.println("temp.size=" + temp.size());
    	curSource = new int[numOfVertice];
    	numOfCurSource = 0;
    	for (int source : temp) {
			curSource[numOfCurSource++] = source;
    	}
    	int c = 50;
    	while (c-- > 0) {
    		getMinFlow();
		}
	}

	private static void verticeInfo() {
		VertexInfo[] capacityInfos = new VertexInfo[numOfVertice];
    	for (int i = 0; i < numOfVertice; i++) {
			int sumCapacity = 0;
			int sumCost = 0;
			for (int index = edgeHeadIndex[i]; index != -1; index = edges[index].next) {
				sumCapacity += edges[index].capacity;
				sumCost += edges[index].capacity * edges[index].cost;
			}
			capacityInfos[i] = new VertexInfo(i, sumCapacity, (double)sumCost / sumCapacity);
		}
    	
    	for (Consumer consumer : consumers) {
			logger.println("id=" + consumer.id + "\tlink=" + consumer.link +
					"\tdemand=" + consumer.demand + "\tvertexInfo=" + capacityInfos[consumer.link]);
		
		}
    	Arrays.sort(capacityInfos, new Comparator<VertexInfo>() {

			@Override
			public int compare(VertexInfo o1, VertexInfo o2) {
				return o1.averageCost < o2.averageCost ? -1 : 1;
			}
		});
    	for (VertexInfo capacityInfo : capacityInfos) {
			logger.println(capacityInfo);
		}
    	logger.println("------------------------------------------------------");
    	logger.flush();
    	
    	//计算贪心算法每种情况的费用
    	for (int num = 1; num <= numOfVertice; num++) {
			int c = 0;
			numOfCurSource = num;
			curSource = new int[numOfCurSource];
			while (c < num) {
				curSource[c] = capacityInfos[c].id;
				c++;
			};
			StringBuilder builder = new StringBuilder(defaultCode);
			for (int i = 0; i < numOfCurSource; i++) {
				builder.setCharAt(curSource[i], '1');
			}
			curCode = builder.toString();
			getMinFlow();
			logger.println("num " + num + " cost: " + curCost);
			LogUtil.printLog("finish num " + num + " cost " + curCost);
		}
	}

	private static void heredity() {
		
		LogUtil.printLog("start heredity");
		
		//参数初始化
//		int numOfPopulation = (((numOfVertice / 2) & 1) == 0) ? 
//				numOfVertice /2 : numOfVertice / 2 +1;
//		int numOfPopulation = ((numOfVertice & 1) == 0) ? numOfVertice : numOfVertice + 1;
//		int numOfPopulation = numOfVertice * 2;
		int numOfPopulation = 50;
		
		String[] population = new String[numOfPopulation];
		int[] fitness = new int[numOfPopulation];
		
		final int stopConst = numOfVertice * 10;
//		final int stopConst = numOfVertice * numOfVertice;
//		final int stopConst = 1000;
		final int stopLine = numOfPopulation / 10;
		int stopCount = 0;
//		final int enhanceLine = 0;
		final int enhanceLine = numOfPopulation * 4 / 10;
		final int maxGenerationConst = 100000;
		
		//消费节点(可以看做和汇点重合)和顶点的有效距离(i为消费节点id，j为顶点下标)
    	int[][] validDistace = new int[numOfConsumers][numOfVertice];
		validDistace(validDistace);
		
		//计算每个顶点之间的距离
//		int[][] distance = new int[numOfVertice][numOfVertice];
//		calDistance(distance);
				
//		int[] maxCapacity = new int[numOfVertice];
//		int maxVertex = -1;
//		int max = -1;
//		for (int vertex = 0; vertex < numOfVertice; vertex++) {
//			int sum = 0;
//			for (int index = edgeHeadIndex[vertex]; index != -1; index = edges[index].next) {
//				sum += edges[index].capacity;
//			}
//			maxCapacity[vertex] = sum;
//			if (sum > max) {
//				maxVertex = vertex;
//				max = sum;
//			}
//		}
		LogUtil.printLog("finish cal validDistance ");
		
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
			
			//层次聚类初始化
			int need = numOfConsumers * 4 / 5;
			int lessThan = numOfConsumers * 4 / 5 + 1;
			List<List<Integer>> situation = levelCluster(validDistace, need, lessThan);
			
			class SourcesCost {
				/**
				 * 源点组合在situation中的id
				 */
				String code;
				int cost;
				public SourcesCost(String code, int cost) {
					this.code = code;
					this.cost = cost;
				}
			}
			
			PriorityQueue<SourcesCost> heap = new PriorityQueue<>(need, new Comparator<SourcesCost>() {

				@Override
				public int compare(SourcesCost o1, SourcesCost o2) {
					return o1.cost - o2.cost;
				}
			});
			
			for (int i = 0; i < need; i++) {
				curSource = new int[numOfVertice];
				numOfCurSource = 0;
				for (int source : situation.get(i)) {
					curSource[numOfCurSource++] = source;
				}
				
				StringBuilder builder = new StringBuilder(defaultCode);
				for (int gene = 0; gene < numOfCurSource; gene++) {
					builder.setCharAt(curSource[gene], '1');
				}
				curCode = builder.toString();
				
				getMinFlow();
				
				if (curCost == Infinity) {
					break;
				}
				
//				LogUtil.printLog("numofsources " + numOfCurSource + " cost " + curCost);
				
				heap.add(new SourcesCost(curCode, curCost));
			}
			
			while (index < numOfPopulation && !heap.isEmpty()) {
				population[index++] = heap.poll().code;
			}
			if (index < numOfPopulation) {
				int len = index;
				int p = 0;
				while (index < numOfPopulation) {
					population[index++] = population[p++];
					p %= len;
				}
			}
			LogUtil.printLog("finish level cluster");
			
			//kmeans聚类初始化种群;
//			int need = numOfPopulation;
//			List<List<Integer>> situation = kMeans(validDistace, need);
//			
//			for (List<Integer> source : situation) {
//				StringBuilder builder = new StringBuilder(defaultCode);
//				for (int gene : source) {
//					builder.setCharAt(gene, '1');
//				}
//				population[index++] = builder.toString();
//			}
//			//防止聚类全部为0,加入极端直连情况
//			StringBuilder builder = new StringBuilder(defaultCode);
//			for (Consumer consumer : consumers) {
//				builder.setCharAt(consumer.link, '1');
//			}
//			population[numOfPopulation - 1] = builder.toString();
//			LogUtil.printLog("finish kmeans");
		
//			index = straightInit(numOfPopulation, population, index);
			
			//随机初始化种群
//			randomInit(numOfPopulation, population, index);

			//计算各个个体的适应度
			fitness = new int[numOfPopulation];
			for (int i = finished; i < numOfPopulation; i++) {
				
//				LogUtil.printLog("cal fintness " + i + "code is " + population[i]);
				
				curCode = population[i];
				if (history.containsKey(curCode)) {
					SourceInfo sourceInfo = history.get(curCode);
					fitness[i] = sourceInfo.fit;
				} else {
					getNewFitness(fitness, i);
				}
			}
			
			for (int fit : fitness) {
				if (fit != 0) {
					inited = true;
					break;
				}
			}
		}		
		
		printOriginalGeneration(population);
	
		LogUtil.printLog("----------------------initTimes: " + initCount);

		//遗传迭代
		int generation = 0;
		for (generation = 0; generation < maxGenerationConst; generation++) {
			if (timer.getTimeDelay() > maxTime) {
				break;
			}
			
//			if (minCost < serverCost * numOfConsumers) {
//				standardCost = serverCost * numOfConsumers / 5 + minCost;
//			}
			
			//计算各个个体的适应度
			fitness = new int[numOfPopulation];
			for (int i = 0; i < numOfPopulation; i++) {
				curCode = population[i];
				if (history.containsKey(curCode)) {
					fitness[i] = history.get(curCode).fit;
				} else {
					getNewFitness(fitness, i);
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

			//保存最优
			
//			int addCount = minCount < 10 ? 10 - minCount : 0;

			int minIndex = 0;
			int minFitness = fitness[0];;
			for (int i = 1; i < fitness.length; i++) {
				if (fitness[i] < minFitness) {
					minIndex = i;
					minFitness = fitness[minIndex];
				}
			}
			//历史最优个体替换当代最差个体
			population[minIndex] = minCode;
			fitness[minIndex] = history.get(minCode).fit;
			minCount++;
		
		
			//适应度尺度变换
			//惩罚非最优个体
			//增大适应度最大的个体的影响力
			int panish = 10;
			int maxFitness = 0;
			double slope= 0;
			
			
//			int panish = minCount < enhanceLine ? (int) (11 - 10.0 / enhanceLine * minCount) : 1;
			if (minCount < enhanceLine) {
				int historyMinFitness = history.get(minCode).fit;
				for (int i = 0; i < population.length; i++) {
					if (fitness[i] != historyMinFitness) {
						fitness[i] /= panish;
					}
				}
				
				minFitness = Infinity;
				for (int fit : fitness) {
					minFitness = Math.min(fit, minFitness);
				}
				maxFitness = history.get(minCode).fit;
				
				if (maxFitness > minFitness) {
					slope = (4.0 * maxFitness - minFitness / 2) / (maxFitness - minFitness);
					for (int i = 0; i < numOfPopulation; i++) {
						fitness[i] = (int) (slope * (fitness[i] - minFitness) + minFitness / 2);
					}
					
				}
			}
			
			
			//输出种群信息
//			if (generation < 3000) {
//				logger.println("max: " + maxFitness + " min: " + minFitness + " slope: " + slope);
//				printpopulation(population, generation, false);
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
//			if (totalFitness < 0) {
//				LogUtil.printLog("totalFitness :" + totalFitness);
//				logger.println("now fitness is :");
//				for (int i =0; i < numOfPopulation; i++) {
//					logger.print(fitness[i] + " ");
//				}
//				logger.println();
//				logger.flush();
//			}
			
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
			
			double crossProbability = 1;
//			double crossProbability = minCount < enhanceLine ? 0.75 : 1;
//			double crossProbability = (0.2 / numOfPopulation) * minCount + 0.8;
			
			//分组交叉
			for (int i = 0; i < numOfPopulation; i += 2) {
//				if (Math.random() > crossProbability) {
//					continue;
//				}
				
				//单点交叉cut,cut+1之间cut[0,len-1);
//				int cut = random.nextInt(nextGeneration[0].length() + 1);
//				String temp = nextGeneration[i].substring(cut);
//				nextGeneration[i] = nextGeneration[i].substring(0, cut) + 
//						nextGeneration[i + 1].substring(cut);
//				nextGeneration[i + 1] = nextGeneration[i + 1].substring(0, cut) + temp;
				
				//两点交叉
				int temp1 = random.nextInt(numOfVertice);
				int temp2 = random.nextInt(numOfVertice);
				int cut1 = Math.min(temp1, temp2);
				int cut2 = Math.max(temp1, temp2);

				String mid = nextGeneration[i].substring(cut1, cut2);
				
				nextGeneration[i] = nextGeneration[i].substring(0, cut1) +
						nextGeneration[i + 1].substring(cut1, cut2) + 
						nextGeneration[i].substring(cut2);
				
				nextGeneration[i + 1] = nextGeneration[i + 1].substring(0, cut1) +
						mid + nextGeneration[i + 1].substring(cut2);
				
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
//				int temp1 = random.nextInt(numOfVertice);
//				int temp2 = random.nextInt(numOfVertice);
//				int cut1 = Math.min(temp1, temp2);
//				int cut2 = Math.max(temp1, temp2);
//				
//				int len = random.nextInt(Math.min(cut2 - cut1, numOfVertice - cut2) + 1);
//				String string = nextGeneration[i];
//				nextGeneration[i] = string.substring(0, cut1) + string.substring(cut2, cut2 + len) +
//						string.substring(cut1 + len, cut2) + string.substring(cut1, cut1 + len) +
//						string.substring(cut2 + len);
//			}
			
			//变异
//			double mutateProbability = minCount < enhanceLine ? 0.1 : 0.2;
			double mutateProbability = (0.15 / numOfPopulation) * minCount + 0.15;
//			double mutateProbability = 0.2;
			
			for (int i= 0; i < numOfPopulation; i++) {

//				char[] genes = nextGeneration[i].toCharArray();
//				if (nextGeneration[i].equals(minCode)) {
//					randomMoveSource(genes);
//				} else if (Math.random() < mutateProbability) {
//					//随机增删源点
//					int pos = random.nextInt(nextGeneration[0].length());
//					genes[pos] = genes[pos] == '0' ? '1' : '0';
//					nextGeneration[i] = new String(genes);
//				}
//				nextGeneration[i] = new String(genes);
				
				if (Math.random() < mutateProbability) {
					char[] genes = nextGeneration[i].toCharArray();
										
					calCost(genes);
					int preCost = curCost;
					String preCode = curCode;
					
					logger.print("before mutate " + preCost + " ; ");

					if (!nextGeneration[i].equals(minCode)) {
//					if (false) {
						Hmove++;
						
						while (true) {
							if (Math.random() < 0.5) {
								//随机移动源点至相邻节点
								randomMoveSource(genes);
							} else {
//								将流量最小的源点移动到相邻的流量最大的顶点
								moveMinFlowToNearlyMaxFlow(genes);
							}
							calCost(genes);
							logger.println("after move " + curCost + " ; change " + (curCost - preCost));
						
							if (curCost >= preCost) {
								genes =  preCode.toCharArray();
							
								break;
							}
							preCost = curCost;
							preCode = curCode;
						}
						
							
					} else if (Math.random() < 0.8) {
						Hdel++;
						
						while (true) {
							if (Math.random() < 0.5) {
								//随机删除源点
								randomDeleteSource(genes);
							} else {
								//删除流量最小的源点
								delMinFlowSource(genes);
							}
							calCost(genes);
							logger.println("after del " + curCost + " ; change " + (curCost - preCost));
							
							if (curCost >= preCost) {
								genes = preCode.toCharArray();
								break;
							}
							preCost = curCost;
							preCode = curCode;
						}
						
					} else {
						Hadd++;
						
						//随机增加源点
						randomAddSource(genes);
						
						calCost(genes);
						logger.println("after add " + curCost + " ; change " + (curCost - preCost));
					}
					
					nextGeneration[i] = new String(genes);
				}
			}
			
			population = nextGeneration;
			
			LogUtil.printLog("generation " + generation + " finish " + 
			"\ncross " + crossProbability + " mutate " + mutateProbability + " min " + minCost);
			
		}
		
		if (logger != null) {
			logger.println("generationCount :" + generation);
		}
		
		LogUtil.printLog("-----------------------generationCount: " + generation + 
				"\nmove " + Hmove + " del " + Hdel + " add " + Hadd);
		
	}

	private static void moveMinFlowToNearlyMaxFlow(char[] genes) {
		int minFlowSource = -1;
		int mf = Infinity;
		for (int j = 0; j < numOfCurSource; j++) {
			int sumflow = 0;
			int source = curSource[j];
//			int sourceOut = source + numOfVertice;
//			for (int k = arcHeadIndex[sourceOut]; k != -1; k = arcs[k].next) {
//				if ((k & 1) == 0) {
//					sumflow += curFlow[k] - curFlow[k ^ 1];
//				}
//			}
			sumflow = curFlow[source * 2] - curFlow[source * 2 + 1];
			
			if (sumflow < mf) {
				minFlowSource = source;
				mf = sumflow;
			}
		}
		
		int maxFlowSource = -1;
		mf = -1;
		for (int j = edgeHeadIndex[minFlowSource]; j != -1; j = edges[j].next) {
			int sumflow = 0;
			int source = edges[j].end;
//			int sourceOut = source + numOfVertice;
//			for (int k = arcHeadIndex[sourceOut]; k != -1; k = arcs[k].next) {
//				if ((k & 1) == 0) {
//					sumflow += curFlow[k] - curFlow[k ^ 1];
//				}
//			}
			sumflow = curFlow[source * 2] - curFlow[source * 2 + 1];
			
			if (sumflow > mf) {
				maxFlowSource = source;
				mf = sumflow;
			}
		}
		
		genes[minFlowSource] = '0';
		genes[maxFlowSource] = '1';
	}

	private static void delMinFlowSource(char[] genes) {
		int ms = -1;
		int mf = Infinity;
		for (int j = 0; j < numOfCurSource; j++) {
			int sumflow = 0;
			int source = curSource[j];
//			int sourceOut = source + numOfVertice;
//			for (int k = arcHeadIndex[sourceOut]; k != -1; k = arcs[k].next) {
//				if ((k & 1) == 0) {
//					sumflow += curFlow[k] - curFlow[k ^ 1];
//				}
//			}
			
			sumflow = curFlow[2 * source] - curFlow[2 * source + 1];
			
			if (sumflow < mf) {
				ms = source;
				mf = sumflow;
			}
		}
		genes[ms] = '0';
		
	}

	private static void calCost(char[] genes) {
		curCode = new String(genes);
		curSource = new int[numOfVertice];
		numOfCurSource = 0;
		for (int gene = 0; gene < numOfVertice; gene++) {
			if (genes[gene] == '1') {
				curSource[numOfCurSource++] = gene;
			}
		}
		int fit = 0;
		if (history.containsKey(curCode)) {
			SourceInfo sourceInfo = history.get(curCode);
			fit = sourceInfo.fit;
			curCost = fit == 0 ? Infinity : (standardCost - fit);
			curFlow = sourceInfo.flow;
		} else {
			getMinFlow();
			fit = curCost > standardCost ? 0 : (standardCost - curCost);
			history.put(curCode, new SourceInfo(fit, curFlow));
		}
		
	}

	private static void randomAddSource(char[] genes) {
		int pos = random.nextInt(numOfVertice);
		while (genes[pos] == '1') {
			pos = random.nextInt(numOfVertice);
		}
		genes[pos] = '1';
	}

	private static void randomDeleteSource(char[] genes) {
		//随机删除源点
		int pos = curSource[random.nextInt(numOfCurSource)];
		genes[pos] = '0';
	}

	private static void randomMoveSource(char[] genes) {
		int from = curSource[random.nextInt(numOfCurSource)];
		int edgeIndex = edgeHeadIndex[from];
		int count = random.nextInt(numOfLinkedEdges[from]);
		while (count-- > 0) {
			edgeIndex = edges[edgeIndex].next;
		}
		int to = edges[edgeIndex].end;
		
		genes[from] = '0';
		genes[to] = '1';
	}
	private static void randomInit(int numOfPopulation, String[] population, int index) {
		while (index < numOfPopulation) {
			int numOfSources = random.nextInt(numOfConsumers) + 1;
			StringBuilder code = new StringBuilder(defaultCode);
			while (numOfSources > 0) {
				int pos = random.nextInt(numOfVertice) ;
				if (code.charAt(pos) == '0') {
					code.setCharAt(pos, '1');
					numOfSources--;
				}
			}
			population[index++] = code.toString();
		}
	}

	private static int straightInit(int numOfPopulation, String[] population, int index) {
		//从极端情况移动源点初始化种群
		curSource = new int[numOfVertice];
		numOfCurSource = 0;
		
		StringBuilder builder2 = new StringBuilder(defaultCode);
		boolean include = true;
		for (Consumer consumer : consumers) {
			if (include) {
				builder2.setCharAt(consumer.link, '1');
				curSource[numOfCurSource++] = consumer.link;
			} 
			include = !include;
		}
		population[index++] = builder2.toString();
		
		while (index < numOfPopulation) {
			int from = curSource[random.nextInt(numOfCurSource)];
			int count = random.nextInt(numOfLinkedEdges[from]);
			int edgeIndex = edgeHeadIndex[from];
			while (count-- > 0) {
				edgeIndex = edges[edgeIndex].next; 
			}
			int to = edges[edgeIndex].end;
			
			builder2.setCharAt(from, '0');
			builder2.setCharAt(to, '1');
			numOfCurSource = 0;
			for (int i = 0; i < builder2.length(); i++) {
				if (builder2.charAt(i) == '1') {
					curSource[numOfCurSource++] = i;
				}
			}
			population[index++] = builder2.toString();
		}
		return index;
	}

	/**
	 * @param population
	 */
	private static void printOriginalGeneration(String[] population) {
		if (logger == null) {
			return;
		}
		List<Integer> sources;
		int validCount = 0;
		logger.println("------------------------");
		logger.println("original generation:");
		int minCost = Infinity;
		for (String unit : population) {
			sources = new ArrayList<>();
			char[] genes = unit.toCharArray();
			for (int gene = 0; gene < genes.length; gene++) {
				if (genes[gene] == '1') {
					sources.add(gene);
				}
			}
			if (history.get(unit).fit != 0) {
				validCount++;
			}
			int cost = history.get(unit).fit == 0 ? -1 : standardCost - history.get(unit).fit;
			logger.println("cost: " + cost +
					" size: " + sources.size() + " " + sources);
			if (cost != -1) {
				minCost = Math.min(minCost, cost);
			}
		}
		logger.println("valid unit amount:" + validCount + " minCost:" + minCost);
		LogUtil.printLog("valid unit amount:" + validCount + " minCost:" + minCost);
		logger.println("-------------------------------");
		logger.flush();
	}

	/**
	 * 遗传算法中获取之前未出现过的个体的适应度
	 * @param fitness
	 * @param position
	 */
	private static void getNewFitness(int[] fitness, int position) {
		numOfCurSource = 0;
		curSource = new int[numOfVertice];
		char[] genes = curCode.toCharArray();
		for (int gene = 0; gene < genes.length; gene++) {
			if (genes[gene] == '1') {
				curSource[numOfCurSource++] = gene;
			}
		}
		
//		for (int i =0; i < numOfCurSource; i++) {
//			System.out.print(curSource[i] + ",");
//		}
//		System.out.println();
		
		getMinFlow();
		
		int fit = curCost > standardCost ? 0 : (standardCost - curCost);
		
		fitness[position] = fit;
		history.put(curCode, new SourceInfo(fit, curFlow));
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
			List<Integer> amount = new ArrayList<>();
			//适应度和amount中的下标
			Map<Integer, Integer> map = new HashMap<>();
			for (String unit : population) {
				int fit = history.get(unit).fit;
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

	private static void sequentialCluster() {
		//消费节点(可以看做和汇点重合)和顶点的有效距离(i为消费节点id，j为顶点下标)
    	int[][] validDistance = new int[numOfVertice][numOfVertice];
		validDistace(validDistance);
		
		int[][] distance = new int[numOfVertice][numOfVertice];
//		calDistance(distance);
		
//		final int maxClusterConst = 100000;
		final int maxClusterConst = 1;
		
//		final int stopConst = 10 * numOfVertice * numOfVertice;
		
		LogUtil.printLog("cluster inited");
		
		int clusterCount = 0;
		int validCount = 0;
		for (clusterCount = 0; clusterCount < maxClusterConst; clusterCount++) {
			if (timer.getTimeDelay() > maxTime) {
				break;
			}
			
			//kmeans获取源点的组合
//			List<List<Integer>> situation = kMeans(validDistance, 1);
			//level获取源点组合
			int need = numOfConsumers * 4 / 5;
			int lessThan = numOfConsumers * 4 / 5 + 1;
			List<List<Integer>> situation = levelCluster(validDistance, need, lessThan);
			
			//计算该源点组合下的最小费用流
	    	for (List<Integer> sources : situation) {
	    		
	    		numOfCurSource = 0;
	    		curSource = new int[numOfVertice];
	    		for (int source : sources) {
					curSource[numOfCurSource++] = source;
				}
	    		
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
				
	    		getMinFlow();
	    		
	    		if (curCost == Infinity) {
					break;
				}
	    		validCount++;
	    		
	    		LogUtil.printLog("numofconsumer " + numOfCurSource + " cost " + curCost);
//		    	if (curCost != Infinity) {
//		    		logger.println(" sourceSize " + numOfCurSource + " cost " + curCost);
//			    }
	    		
	    		//记录计算过最小费用流的编码
	    		history.put(curCode, new SourceInfo(0, curFlow));
	    	}
	    	
	    }
		
		LogUtil.printLog("--------------clusterCount: " + clusterCount + " validCount " + validCount);
	}

	private static void calDistance(int[][] distance) {
		for (int v0 = 0; v0 < numOfVertice; v0++) {
			Path[] paths = dijstra(v0);
			for (int j = 0; j < numOfVertice; j++) {
				distance[v0][j] = paths[j].distance;
			}
		}
	}
    
	/**
	 * 通过最小费用流矩阵求得路径
	 * @param graphO
	 * @param edges
	 */
	private static void getResult() {
		//增加超级源点,超级源点和源点的弧及对应的反向弧
//		int tempNumOfArcs = numOfArcs;
		
		for (int i = 0; i < numOfMinSource; i++) {
			int sourceIn = minSource[i];
			//<out1,in2>及其反向弧
			arcs[numOfArcs] = new Arc(superSource, sourceIn, Infinity, 0, arcHeadIndex[superSource]);
			arcHeadIndex[superSource] = numOfArcs++;
//			numOfLinkedArcs[superSource]++;
			//arcIndex奇数表示为反向弧
			arcs[numOfArcs] = new Arc(sourceIn, superSource, 0, 0, arcHeadIndex[sourceIn]);
			arcHeadIndex[sourceIn] = numOfArcs++;
//			numOfLinkedArcs[sourceOut]++;
		}
    	
		//用反向弧的流量撤销正向弧的流量
		for (int i = 0; i < numOfArcs; i+= 2) {
			minFlow[i] -= minFlow[i + 1];
			minFlow[i + 1] = 0;
		}
		
		//合并出点和入点，去除反向弧，转化流网络
		//或者说，在原始网络上增加超级源汇点和对应弧
		Arc[] tempEdges = edges;
		int[] tempedgeHeadIndex = edgeHeadIndex;
		
		edges = new Arc[2 * numOfEdges + numOfMinSource + numOfConsumers];
		edgeHeadIndex = new int[numOfVertice + 2];
		
		System.arraycopy(tempEdges, 0, edges, 0, tempEdges.length);
		System.arraycopy(tempedgeHeadIndex, 0, edgeHeadIndex, 0, tempedgeHeadIndex.length);
		edgeHeadIndex[numOfVertice] = -1;
		edgeHeadIndex[numOfVertice + 1] = -1;
		
		//超级源点
		int edgeIndex = tempEdges.length;
		for (int i = 0; i < numOfMinSource; i++) {
			edges[edgeIndex] = new Arc(numOfVertice, minSource[i], Infinity, 0, edgeHeadIndex[numOfVertice]);
			edgeHeadIndex[numOfVertice] = edgeIndex++;
		}
		//超级汇点
		for (Consumer consumer : consumers) {
			int sink = consumer.link;
			edges[edgeIndex] = new Arc(sink, numOfVertice + 1, consumer.demand, 0, edgeHeadIndex[sink]);
			edgeHeadIndex[sink] = edgeIndex++;
		}
		
//		for (int i = 0; i < numOfArcs; i++) {
//			if (minFlow[i] != 0) {
//				logger.println(i + " " + minFlow[i]);
//			}
//		}
		
		int[][] remain = new int[numOfVertice + 2][numOfVertice + 2];
		
		for (int i = 2 * numOfVertice; i < numOfArcs; i += 2) {
			int start = arcs[i].start;
			int end = arcs[i].end;
			if (end == superSink) {
				end = numOfVertice + 1;
			}
			int flow = minFlow[i];
			remain[start - numOfVertice][end] = flow;
		}
		
		//保证超级源点的流量充足
		for (int i = 0; i < numOfMinSource; i++) {
			remain[numOfVertice][minSource[i]] = Infinity;
		}
		
		
//		for (int i = 0; i < remain.length; i++) {
//			for (int j = 0; j < remain[0].length; j++) {
//				if (remain[i][j] == 0) {
//					continue;
//				}
//				if (remain[i][j] == remain[j][i]) {
//					logger.println("[" + i + "][" + j + "]=" + remain[i][j] + ",");
//				}
//			}
//		}
		
//		printRemain(remain);
		
//		for (int i = 0; i < numOfMinSource; i++) {
//			int source = minSource[i];
//			int sumFlow = 0;
//			int sumCapacity = 0;
//			logger.println("source " + source + " :");
//			for (int k = edgeHeadIndex[source]; k != -1; k = edges[k].next) {
//				int flow = remain[source][edges[k].end];
//				int capacity = edges[k].capacity;
//				sumFlow += flow;
//				sumCapacity += capacity;
//				logger.println("flow " + flow + " capacity " + capacity + " cost " + edges[k].cost);
//				
//			}
//			logger.println("sumFlow " + sumFlow + " sumcapacity " + sumCapacity);
//		}
//		logger.flush();

		
    	while (true) {
			if (!dfs(new boolean[numOfVertice + 2], numOfVertice, Infinity, new ArrayList<Integer>(numOfConsumers), remain)) {
				break;
			}
		}
    	
    	int count = result.size();
    	
    	
//		printRemain(remain);

    	result.addAll(0, Arrays.asList(new String[]{Integer.toString(count),""}));
	}

	private static void printRemain(int[][] remain) {
		logger.println("remain flow");
		for (int i = 0; i < numOfVertice + 2; i++) {
			for (int j = 0; j < numOfVertice + 2; j++) {
				logger.print(remain[i][j] + " ");
			}
			logger.println();
		}
		logger.println("remain[26][0]=" + remain[26][0]);
		logger.println("------------------------------------------");
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
			}
			return true;
		}
    	
    	path.add(cur);
    	visited[cur] = true;
    	
    	for (int edgeIndex = edgeHeadIndex[cur]; edgeIndex != -1; edgeIndex = edges[edgeIndex].next) {
			int next = edges[edgeIndex].end;
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
	/**
	 * 根据源点组合建立有向图，求最小费用流
	 * @param graphO
	 * @param result
	 * @param sources
	 */
	private static void getMinFlow() {
		
		long time1 = System.currentTimeMillis();
//		LogUtil.printLog("start minFlow");
		
		//增加超级源点,超级源点和源点的弧及对应的反向弧
		int tempNumOfArcs = numOfArcs;
		
		for (int i = 0; i < numOfCurSource; i++) {
			int sourceIn = curSource[i];
			//<out1,in2>及其反向弧
			arcs[numOfArcs] = new Arc(superSource, sourceIn, Infinity, 0, arcHeadIndex[superSource]);
			arcHeadIndex[superSource] = numOfArcs++;
//			numOfLinkedArcs[superSource]++;
			//arcIndex奇数表示为反向弧
			arcs[numOfArcs] = new Arc(sourceIn, superSource, 0, 0, arcHeadIndex[sourceIn]);
			arcHeadIndex[sourceIn] = numOfArcs++;
//			numOfLinkedArcs[sourceOut]++;
		}
//		LogUtil.printLog("finish graphDi");
		
		//有向图的流网络
		int[] flowDi = new int[numOfArcs];
		//有向图的剩余网络
		int[]remain = new int[numOfArcs];
		
		//初始化流网络
		for (int i = 0; i < numOfArcs; i+= 2) {
			remain[i] = arcs[i].capacity;
		}
		
//		LogUtil.printLog("start expand");

		int queueSize = numOfVertice;
//		int[] arcHeadIndex = Deploy.arcHeadIndex;
//		Arc[] arcs = Deploy.arcs;
		int count = 0;
		
		int[] preArc = null;
		int[] distance = null;
		boolean done = false;
		
		while (true) {
			//SPFA计算最短路径
//			System.out.println(++count);
			//前驱节点下标
//			int[] preVertex = new int[diNumOfVertice];
//			Arrays.fill(preVertex, -1);
			//前驱弧在arcs中的下标
			preArc = new int[diNumOfVertice];
			Arrays.fill(preArc, -1);
			//到superSource的距离
			distance = new int[diNumOfVertice];
			Arrays.fill(distance, Infinity);
			
			distance[superSource] = 0;
			
			//节点是否在队列中
			boolean[] inQueue = new boolean[diNumOfVertice];
			//循环队列效果比大空间的队列效果好
			int[] queue = new int[queueSize]; 
			int head = 0;
			int tail = 0;
			queue[tail++] = superSource;
			inQueue[superSource] = true;
			
//			int[] outQueueCount = new int[diNumOfVertice];
//			int sumDistaceInqueue = 0;
			
//			int queueCount = 0;
			
			while (head != tail) {
//				if (count == 1257) {
//					System.out.println("queueSize= " + (tail + queueSize - head) % queueSize);
//				}
				
				//LLL优化，慢5%-10%????
//				while (distance[queue[head]] * (tail - head) > sumDistaceInqueue) {
//					queue[tail++] = queue[head++];
//				}
				int start = queue[head++];
				if (head == queueSize) {
					head = 0;
				}
				inQueue[start] = false;
				
//				if (++outQueueCount[start] > 10) {
//					done = true;
//					break;
//				}
//				System.out.println("outqueueCount " + start + " = " + outQueueCount[start]);
				
				for (int arcIndex = arcHeadIndex[start]; arcIndex != -1; arcIndex = arcs[arcIndex].next) {
					if (remain[arcIndex] > 0) {
						int end = arcs[arcIndex].end;
						int cost= arcs[arcIndex].cost;
						
						if (distance[start] + cost < distance[end]) {
							distance[end] = distance[start] + cost;
//							preVertex[end] = start;
							preArc[end] = arcIndex;
							if (!inQueue[end]) {
								//SLF优化
								if (distance[end] < distance[queue[head]]) {
									if (head == 0) {
										head = queueSize;
									}
									queue[--head] = end;
								} else {
									queue[tail++] = end;
									if (tail == queueSize) {
										tail = 0;
									}
								}
								inQueue[end] = true;
//								sumDistaceInqueue += distance[end];
							}
						}
					}
				}
			}
			
			//无最短路径，已获得最小费用流
			if (distance[superSink] == Infinity || done) {
				break;
			}
			
			//最短路径的流量
			int minFlow = Infinity;
			for (int arcIndex = preArc[superSink]; arcIndex != -1; arcIndex = preArc[arcs[arcIndex].start] ) {
				minFlow = Math.min(minFlow, remain[arcIndex]);
			}
			
			//求剩余网络
			for (int arcIndex = preArc[superSink]; arcIndex != -1; arcIndex = preArc[arcs[arcIndex].start]) {
				flowDi[arcIndex] += minFlow;
				remain[arcIndex] -= minFlow;
				remain[arcIndex ^ 1] += minFlow;
			}
			
		}
		
//		LogUtil.printLog("finish expand");
		
		//是否满足消费节点的需求
		int sum = 0;
		for (int arcIndex = arcHeadIndex[superSource]; arcIndex != -1; arcIndex = arcs[arcIndex].next) {
			sum += flowDi[arcIndex];
		}
		boolean enough = sum == totalDemand;
		
		int sumCost = Infinity;
		if (enough) {
			sumCost = 0;
			
			for (int i = 0; i < numOfArcs; i += 2) {
				sumCost += (flowDi[i] - flowDi[i + 1]) * arcs[i].cost;
			}
			
			sumCost += serverCost * numOfCurSource;
			
			if (sumCost < minCost) {
				minCost = sumCost;
				minFlow = flowDi;
				minSource = curSource;
				numOfMinSource = numOfCurSource;
				minCode = curCode;
				
//				boolean sourceFull = true;
//				for (int k = arcHeadIndex[superSource]; k != -1; k = arcs[k].next) {
//					if (minFlow[k] == 0) {
//						sourceFull = false;
//					}
//				}
//				logger.println("minFlow full : " + sourceFull);
			}
			

		}

		curCost = sumCost;
		curFlow = flowDi;
		
//		LogUtil.printLog("prepare for next");
		
		//删除有向图超级源点到源点的弧及其反向弧
		for (int i = 0; i < numOfCurSource; i++) {
			int sourceIn = curSource[i];
			arcHeadIndex[sourceIn] = arcs[arcHeadIndex[sourceIn]].next;
		}
		arcHeadIndex[superSource] = -1;
		
		numOfArcs = tempNumOfArcs;
		
//		System.out.println(System.currentTimeMillis() - time1);
		getMinTime += System.currentTimeMillis() - time1;
		getMinCount++;
//		LogUtil.printLog("finish minFlow");
		
	}

	/**
	 * k-means++聚类算法
	 * @param distance
	 * @param need
	 * @return
	 */
//	private static List<List<Integer>> kMeansPP(int[][] distance, int need) {
//		
//	}
	
	/**
	 * k-means聚类寻找源点组合
	 * @param distance
	 * @param need 需要的源点组合个数
	 * @return
	 */
	private static List<List<Integer>> kMeans(int[][] validDistance, int need) {
		//源点组合
    	List<List<Integer>> situation = new ArrayList<>(numOfConsumers);
    	
    	//开始聚类时存在num个源点
    	while (need-- > 0) {
			int num = random.nextInt(numOfConsumers * 2 / 5) + numOfConsumers * 2 / 5;
//    		int num = numOfConsumers / 2;
    		//源点的位置
			List<Integer> newSources = new ArrayList<>(num);
			
			//随机获得num个源点
			while (newSources.size() < num) {
				int temp = random.nextInt(numOfVertice);
				if (!newSources.contains(temp)) {
					newSources.add(temp);
				}
			}
			
			//随机初始源点，之后选取距离已有节点最小距离最大的节点（汇点作为源点）
			//记录每个节点距离已有节点的最小距离
//			int[] closest = new int[numOfVertice];
//			//消费节点的下标
//			int first = random.nextInt(numOfVertice);
//			newSources.add(first);
//			for (int i = 0; i < numOfVertice; i++) {
//				closest[i] = distace[first][i];
//			}
//			closest[first] = -1;
//			
//			while (newSources.size() < num) {
//				int temp = 0;
//				int max = closest[temp];
//				for (int i = 1; i < numOfVertice; i++) {
//					if (closest[i] > max) {
//						temp = i;
//						max = closest[temp];
//					}
//				}
//				//选取的节点季倩如初始源点组合
//				newSources.add(temp);
//				//更新closest
//				for (int i = 0; i < numOfVertice; i++) {
//					closest[i] = Math.min(closest[i], distace[temp][i]);
//				}
//				closest[temp] = -1;
//			}
			
			//K=源点下标 V=属于该源点的消费节点id的集合
			Map<Integer, List<Integer>> map = new HashMap<>(num);
			
			//寻找距离汇点最近的源点，则该汇点属于该源点
			for (int consumerId = 0; consumerId < numOfConsumers; consumerId++) {
				int minSource = -1;
				int minDistance = Infinity;
				for (int source : newSources) {
					if (validDistance[consumerId][source] < minDistance) {
						minSource = source;
						minDistance = validDistance[consumerId][minSource];
					}
				}
				if (!map.containsKey(minSource)) {
					map.put(minSource, new ArrayList<Integer>(numOfConsumers));
				}
				map.get(minSource).add(consumerId);
			}
			
			boolean move = true;
			
			//标准k-means，效果不好，原因未知。。。
			while (move) {

//				if (logger != null) {
//					logger.println(map.keySet());
//				}
				
				move = false;
				
				//每一组汇点寻找新的中心作为新的源点
				newSources = new ArrayList<>(num);
				for (int source : map.keySet()) {
					int newCenter = -1;
					int minSum = Infinity;
					for (int vertex = 0; vertex < numOfVertice; vertex++) {
						int sum = 0;
						for (int consumerId : map.get(source)) {
							sum += validDistance[consumerId][vertex];
						}
						if (sum < minSum) {
							newCenter = vertex;
							minSum = sum;
						}
					}
					
					if (!newSources.contains(newCenter)) {
						newSources.add(newCenter);
					}
					
					if (source != newCenter) {
						move = true;
					} 
				}
				
				map = new HashMap<>(num);
				
				//寻找距离汇点最近的源点，则该汇点属于该源点
				for (int consumerId = 0; consumerId < numOfConsumers; consumerId++) {
					int minSource = -1;
					int minDistance = Infinity;
					for (int source : newSources) {
						if (validDistance[consumerId][source] < minDistance) {
							minSource = source;
							minDistance = validDistance[consumerId][minSource];
						}
					}
					if (!map.containsKey(minSource)) {
						map.put(minSource, new ArrayList<Integer>(numOfConsumers));
					}
					map.get(minSource).add(consumerId);
				}
			}
			
//			if (logger != null) {
//				logger.println(map.keySet());
//				logger.println("----------------------------------------");
//			}
			
			
			//改编的k-means
			//初始化在一个集合的汇点直到最终也在同一个集合
//			while (move) {
//				move = false;
//				Map<Integer, List<Integer>> newMap = new HashMap<>(num);
//				
//				//对每一组汇点选择新的中心作为新的源点
//				for (int source : map.keySet()) {
//					List<Integer> sinkCollection = map.get(source);
//					int newCenter = -1;
//					int minSum = Infinity;
//					for (int vertex = 0; vertex < numOfVertice; vertex++) {
//						int sum = 0;
//						for (int consumerId : sinkCollection) {
//							sum += validDistance[consumerId][vertex];
//						}
//						if (sum < minSum) {
//							newCenter = vertex;
//							minSum = sum;
//						}
//					}
//					if (source != newCenter) {
//						move = true;
//					} 
//					
//					//判断是否与其他组的源点相同
//					if (newMap.containsKey(newCenter)) {
//						newMap.get(newCenter).addAll(sinkCollection);
//					} else {
//						newMap.put(newCenter, sinkCollection);
//					}
//				}
//				map = newMap;
//			}
			
			situation.add(new ArrayList<>(map.keySet()));
		
		}
    	
    	return situation;
	}

	/**
	 * 原始网络中各个顶点到汇点的距离
	 * @param edges
	 * @param consumerDistance
	 */
	private static void validDistace(int[][] consumerDistance) {
		int minDemand = Infinity;
		for (Consumer consumer : consumers) {
			minDemand = Math.min(minDemand, consumer.demand);
		}
		
		for (int id = 0; id < numOfConsumers; id++) {
    		Consumer consumer = consumers[id];
			int sink = consumer.link;
			
			Path[] paths = dijstra(sink);
			
			for (int i = 0; i < numOfVertice; i++) {
//				distance[id][i] = (paths[i].distance << 16) / consumer.demand;
				
				int weight = (consumer.demand - minDemand) * 10 + minDemand;
				
				consumerDistance[id][i] = paths[i].distance * weight;
			}
			
		}
		
//		if (logger != null) {
//			logger.println("shortestDistance:");
//			logger.println("----------------------------------------------------");
//			for (int i = 0; i < consumerDistance.length; i++) {
//				for (int j = 0; j < consumerDistance[0].length; j++) {
//					logger.print(consumerDistance[i][j] + " ");
//				}
//				logger.println();
//			}
//			logger.println("----------------------------------------------------");
//		}
		
	}

	/**
	 * heap + dijstra
	 * @param edges
	 * @param v0 最短路径的起点
	 * @return
	 */
	private static Path[] dijstra(int v0) {
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
				paths[vertex] = new Path(vertex, -1, getEdgeCost(v0, vertex));
//				paths[vertex] = new Path(vertex, -1, Infinity);
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
			for (int edgeIndex = edgeHeadIndex[start]; edgeIndex != -1; edgeIndex = edges[edgeIndex].next) {
				Arc edge = edges[edgeIndex];
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

	private static int getEdgeCost(int start, int end) {
		for (int edgeIndex = edgeHeadIndex[start]; edgeIndex != -1; edgeIndex = edges[edgeIndex].next) {
			Arc edge = edges[edgeIndex];
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
    	numOfVertice = Integer.parseInt(tokens[0]);
    	numOfEdges = Integer.parseInt(tokens[1]);
    	numOfConsumers = Integer.parseInt(tokens[2]);
    	serverCost = Integer.parseInt(graphContent[2]);
    	
    	consumers = new Consumer[numOfConsumers];
    	getConsumerId = new HashMap<>(numOfConsumers);
    	diNumOfVertice = numOfVertice * 2 + 2;
    	superSource = numOfVertice * 2;
    	superSink = superSource + 1;
    	
    	edges = new Arc[2 * numOfEdges];
    	edgeHeadIndex = new int[numOfVertice];
    	Arrays.fill(edgeHeadIndex, -1);
    	numOfLinkedEdges = new int[numOfVertice];
    	
    	arcs = new Arc[ numOfVertice * 2 + numOfEdges * 4 + 
    	                numOfConsumers * 2 + numOfVertice * 2];
    	arcHeadIndex = new int[diNumOfVertice];
    	Arrays.fill(arcHeadIndex, -1);
//    	numOfLinkedArcs = new int[diNumOfVertice];
    	
    	defaultCode = new StringBuilder(numOfVertice);
    	for (int i = 0; i < numOfVertice; i++) {
			defaultCode.append('0');
		}
    	
    	history = new HashMap<>(numOfVertice * numOfVertice);
    	history.put(defaultCode.toString(), new SourceInfo(0, null));
    	
    	standardCost = serverCost * numOfConsumers * 2;
    	
    	//初始化有向图的弧和对应的反向弧，除了超级源点和源点的弧
    	
    	int dataIndex = 0;
    	int arcIndex = 0;
    	int edgeIndex= 0;
    	
    	//节点内部的弧及其反向弧
    	for (int in = 0; in < numOfVertice; in++) {
			int out = in + numOfVertice;
			
			//内部弧
			arcs[arcIndex] = new Arc(in, out, Infinity, 0, arcHeadIndex[in]);
			arcHeadIndex[in] = arcIndex++;
//			numOfLinkedArcs[in]++;
			
			//反向弧
			arcs[arcIndex] = new Arc(out, in, 0, 0, arcHeadIndex[out]);
			arcHeadIndex[out] = arcIndex++;
//			numOfLinkedArcs[out]++;
		}
    	
//    	logger.println("split vertex");
//    	logger.println("arcNum " + arcIndex + " edgeNum " + edgeIndex);
    	
    	for (dataIndex = 4; dataIndex < graphContent.length; dataIndex++) {
    		String line = graphContent[dataIndex];
    		if (line.length()==0) {
				break;
			}
			String[] edge = line.split(" ");
			int v1 = Integer.parseInt(edge[0]);
			int v2 = Integer.parseInt(edge[1]);
			int capacity = Integer.parseInt(edge[2]);
			int cost = Integer.parseInt(edge[3]);
			
			edges[edgeIndex] = new Arc(v1, v2, capacity, cost, edgeHeadIndex[v1]);
			edgeHeadIndex[v1] = edgeIndex++;
			numOfLinkedEdges[v1]++;
			edges[edgeIndex] = new Arc(v2, v1, capacity, cost, edgeHeadIndex[v2]);
			edgeHeadIndex[v2] = edgeIndex++;
			numOfLinkedEdges[v2]++;
			
			int in1 = v1;
			int out1 = in1 + numOfVertice;
			int in2 = v2;
			int out2 = in2 + numOfVertice;
			
			//<out1,in2>及其反向弧
			arcs[arcIndex] = new Arc(out1, in2, capacity, cost, arcHeadIndex[out1]);
			arcHeadIndex[out1] = arcIndex++;
//			numOfLinkedArcs[out1]++;
			//arcIndex奇数表示为反向弧
			arcs[arcIndex] = new Arc(in2, out1, 0, -cost, arcHeadIndex[in2]);
			arcHeadIndex[in2] = arcIndex++;
//			numOfLinkedArcs[in2]++;
			
			//<out2,in1>及其反向弧
			arcs[arcIndex] = new Arc(out2, in1, capacity, cost, arcHeadIndex[out2]);
			arcHeadIndex[out2] = arcIndex++;
//			numOfLinkedArcs[out2]++;
			arcs[arcIndex] = new Arc(in1, out2, 0, -cost, arcHeadIndex[in1]);
			arcHeadIndex[in1] = arcIndex++;
//			numOfLinkedArcs[in1]++;
		}
//    	
//    	logger.println("split edge");
//    	logger.println("arcNum " + arcIndex + " edgeNum " + edgeIndex);
    	
    	//初始化消费节点和连接的节点（汇点），汇点到超级汇点的弧及其反向弧
    	int conIndex = 0;
    	for (dataIndex++; dataIndex < graphContent.length; dataIndex++) {
			String[] consumer = graphContent[dataIndex].split(" ");
			int id = Integer.parseInt(consumer[0]);
			int link = Integer.parseInt(consumer[1]);
			int demand = Integer.parseInt(consumer[2]);
			
			totalDemand += demand;
			
			int sinkOut = link + numOfVertice;
			consumers[conIndex++] = new Consumer(id, link, demand);
			getConsumerId.put(link, id);
			
			arcs[arcIndex] = new Arc(sinkOut, superSink, demand, 0, arcHeadIndex[sinkOut]);
			arcHeadIndex[sinkOut] = arcIndex++;
//			numOfLinkedArcs[sinkOut]++;
			arcs[arcIndex] = new Arc(superSink, sinkOut, 0, 0, arcHeadIndex[superSink]);
			arcHeadIndex[superSink] = arcIndex++;
//			numOfLinkedArcs[superSink]++;
		}
    	

//    	logger.println("add consumers");
//    	logger.println("arcNum " + arcIndex + " edgeNum " + edgeIndex);
//    	logger.flush();
    	
		numOfArcs = arcIndex;
		
	}
    
	private static void printArcAndEdge(boolean detail) {
		int arcIndex;
		int edgeIndex;
		logger.println("edge info");
		logger.println("numofedges " + numOfEdges);
		if (detail) {
			for (int vertex = 0; vertex < numOfVertice; vertex++) {
				for (edgeIndex = edgeHeadIndex[vertex]; edgeIndex != -1; edgeIndex = edges[edgeIndex].next) {
					Arc edge = edges[edgeIndex];
					logger.println(edge);
				}
			}
		}
		logger.println("-----------------------------------------");
		logger.println("numofarcs " + numOfArcs);
		logger.println("arc info");
		if (detail) {
			for (int vertex = 0; vertex <= superSink; vertex++) {
				for (arcIndex = arcHeadIndex[vertex]; arcIndex != -1; arcIndex = arcs[arcIndex].next) {
					Arc arc = arcs[arcIndex];
					logger.println(arc);
				}
			}
		}
		logger.println("-----------------------------------------");
		logger.flush();
	}
    
    
	private static class SourceInfo {
		private int fit;
		private int[] flow;
		
		public SourceInfo(int fit, int[] flow) {
			this.fit = fit;
			this.flow = flow;
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

	/**
	 * 有向图中的弧，是否为反向弧通过数组中的位置是否为奇数判断
	 * @author hjg
	 *
	 */
	private static class Arc {
		
		private int start;
		
		private int end;
		
		private int capacity;
		
		private int cost;
		
		/**
		 * 下一个节点在arcs中的下标,-1表示无下一个节点
		 */
		private int next;
		
		public Arc(int start, int end, int capacity, int cost, int next) {
			this.start = start;
			this.end = end;
			this.capacity = capacity;
			this.cost = cost;
			this.next = next;
		}
		
		@Override
		public String toString() {
			return "["+start+"->"+end+": capacity "+capacity+" cost "+cost + "]";
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
	
	public static class VertexInfo {
		private int id;
		private int sumFlow;
		private double averageCost;
		
		public VertexInfo(int id, int sumFlow, double averageCost) {
			this.id  = id;
			this.sumFlow = sumFlow;
			this.averageCost = averageCost;
		}
		
		@Override
		public String toString() {
			return "[id=" + id +",sum=" + sumFlow + ",averageCost=" + averageCost + "]";
		}
	}
}
