import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.io.*;

class TSP {
	static final int MAX_N = 1005;

	static int N, position, bestTourDistance;
	static long end;
	static boolean hasCycle;
	static int[] tour, bestTour, cycleParent;
	static boolean[] visited;
	static double[][] points = new double[MAX_N][2];
	static int[][] distance = new int[MAX_N][MAX_N], memo, parent;
	static ArrayList<ArrayList<Integer>> adjList;
	static PriorityQueue<IntegerTriple> edgeList;

	public static void main(String[] args) throws Exception {
		MyScanner sc = new MyScanner();
		PrintWriter out = new PrintWriter(new BufferedOutputStream(System.out));

		while (true) {
			String line = sc.nextLine();
			if (line == null) {
				break;
			}

			long start = System.currentTimeMillis();
			end = start + 1550;

			N = Integer.parseInt(line);
			tour = new int[N];
			bestTourDistance = Integer.MAX_VALUE;
			edgeList = new PriorityQueue<IntegerTriple>();

			for (int i = 0; i < N; i++) {
				points[i][0] = sc.nextDouble();
				points[i][1] = sc.nextDouble();
			}

			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					distance[j][i] = distance[i][j] = EuclideanDistance(points[j][0] - points[i][0],
							points[j][1] - points[i][1]);
					edgeList.add(new IntegerTriple(distance[i][j], i, j));
				}
			}

			if (N < 1) {
				dp();
			} else {
				MultipleFragmentTour();

				while (System.currentTimeMillis() < end) {
					twoOpt();
					twoHalfOpt();
					CompareTour();
					LocalRandom();
					// DoubleBridge();
				}
			}

			CompareTour();
			for (int i = 0; i < N; i++) {
				out.println(bestTour[i]);
			}
		}

		out.close();
	}

	static int EuclideanDistance(double x, double y) {
		double distance = Math.sqrt((x * x) + (y * y));
		return (int) Math.round(distance);
	}

	static void dp() {
		memo = new int[N][1 << N];
		parent = new int[N][1 << N];
		Query(0, 1);

		int path_counter = 0;
		int cur_node = 0;
		int cur_mask = 1;

		do {
			tour[path_counter] = cur_node;
			path_counter++;
			cur_node = parent[cur_node][cur_mask];
			cur_mask = cur_mask | (1 << cur_node);
		} while (cur_node - 1 != -1);
	}

	static int Query(int u, int m) {
		if (m == (1 << N) - 1)
			return distance[u][0];

		if (memo[u][m] - 1 != -1)
			return memo[u][m];

		memo[u][m] = Integer.MAX_VALUE;
		for (int i = 0; i < N; i++) {
			if (i != u && (m & (1 << i)) == 0) {
				int cur_value = distance[u][i] + Query(i, m | (1 << i));
				if (memo[u][m] > cur_value) {
					memo[u][m] = cur_value;
					parent[u][m] = i;
				}
			}
		}

		return memo[u][m];
	}

	static void MultipleFragmentTour() {
		adjList = new ArrayList<ArrayList<Integer>>();

		for (int i = 0; i < N; i++) {
			adjList.add(new ArrayList<Integer>());
		}

		int[] taken = new int[N];
		visited = new boolean[N];
		cycleParent = new int[N];
		Arrays.fill(cycleParent, -1);

		while (!edgeList.isEmpty()) {
			IntegerTriple top = edgeList.poll();
			int v = top.second();
			int u = top.third();

			if (taken[v] != 2 && taken[u] != 2) {
				adjList.get(v).add(u);
				adjList.get(u).add(v);

				Arrays.fill(visited, false);
				hasCycle = false;
				CheckCycle(v);

				if (hasCycle) {
					adjList.get(v).remove((Integer) u);
					adjList.get(u).remove((Integer) v);
				} else {
					taken[v]++;
					taken[u]++;
				}
			}
		}

		Arrays.fill(visited, false);
		position = 0;
		DFS(0);
	}

	static void CheckCycle(int v) {
		visited[v] = true;

		ArrayList<Integer> neighbors = adjList.get(v);
		for (int i = 0; i < neighbors.size(); i++) {
			int u = neighbors.get(i);
			if ((visited[u] && cycleParent[v] != u) || hasCycle) {
				hasCycle = true;
				return;
			} else if (!visited[u]) {
				cycleParent[u] = v;
				CheckCycle(u);
			}
		}
	}

	static void DFS(int v) {
		visited[v] = true;
		tour[position++] = v;

		ArrayList<Integer> neighbors = adjList.get(v);
		for (int i = 0; i < neighbors.size(); i++) {
			int u = neighbors.get(i);
			if (!visited[u]) {
				DFS(u);
			}
		}
	}

	static void twoOpt() {
		int minChange, a, b, c, d;

		do {
			if (System.currentTimeMillis() > end)
				return;

			minChange = 0;
			int mini = -1, minj = -1, change;

			for (int i = 0; i < N - 2; i++) {
				a = tour[i];
				b = tour[i + 1];

				for (int j = i + 2; j < N; j++) {
					c = tour[j];
					d = (j == N - 1) ? tour[0] : tour[j + 1];
					change = distance[a][c] + distance[b][d] - distance[a][b] - distance[c][d];

					if (minChange > change) {
						minChange = change;
						mini = i;
						minj = j;
					}
				}
			}

			if (mini != -1 && minj != -1)
				twoOptSwap(mini, minj);

		} while (minChange < 0);
	}

	static void twoOptSwap(int first, int second) {
		int limit, pos, pos2;
		if (second - first <= N / 2) {
			limit = (second - first) / 2;
			pos = first + 1;
			pos2 = second;
		} else {
			limit = ((N - second - 1) + (first + 1)) / 2;
			pos = (second == N - 1) ? 0 : second + 1;
			pos2 = first;
		}

		for (int i = 0; i < limit; i++) {
			int tmp = tour[pos];
			tour[pos++] = tour[pos2];
			tour[pos2--] = tmp;

			if (pos == N)
				pos = 0;
			if (pos2 == -1)
				pos2 = N - 1;
		}
	}

	static void twoHalfOpt() {
		int minChange, a, b, c, d, e;

		do {
			if (System.currentTimeMillis() > end)
				return;

			minChange = 0;
			int mini = -1, minj = -1, change;

			for (int i = 0; i < N - 3; i++) {
				a = tour[i];
				b = tour[i + 1];
				c = tour[i + 2];

				for (int j = i + 3; j < N; j++) {
					d = tour[j];
					e = (j == N - 1) ? tour[0] : tour[j + 1];
					change = distance[d][b] + distance[b][e] + distance[a][c] - distance[a][b] - distance[b][c]
							- distance[d][e];

					if (minChange > change) {
						minChange = change;
						mini = i;
						minj = j;
					}
				}
			}

			if (mini != -1 && minj != -1)
				twoHalfOptSwap(mini, minj);

		} while (minChange < 0);
	}

	static void twoHalfOptSwap(int first, int second) {
		int tmp = tour[first + 1];

		for (int i = first + 1; i < second; i++) {
			tour[i] = tour[i + 1];
		}

		tour[second] = tmp;
	}

	static void CompareTour() {
		int distance = totalDistance();
		if (distance < bestTourDistance) {
			bestTourDistance = distance;
			bestTour = tour.clone();
		} else {
			tour = bestTour.clone();
		}
	}

	static int totalDistance() {
		int total = distance[tour[N - 1]][tour[0]];

		for (int i = 0; i < N - 1; i++) {
			total += distance[tour[i]][tour[i + 1]];
		}

		return total;
	}

	static void LocalRandom() {
		int K = N / 15;
		for (int i = 0; i < K; i++) {
			int pos = ThreadLocalRandom.current().nextInt(N);
			int pos2 = (pos ^ (pos ^ K)) % N;
			int tmp = tour[pos];
			tour[pos] = tour[pos2];
			tour[pos2] = tmp;
		}
	}

	static void DoubleBridge() {
		int random = ThreadLocalRandom.current().nextInt(N / 4) + 1;
		int A = ThreadLocalRandom.current().nextInt(N / 4) + 1;
		int B = A + random;
		int C = B + random;

		int pos = 0;
		int[] tmp = tour.clone();
		for (int i = 0; i < A; i++)
			tour[pos++] = tmp[i];

		for (int i = C; i < N; i++)
			tour[pos++] = tmp[i];

		for (int i = B; i < C; i++)
			tour[pos++] = tmp[i];

		for (int i = A; i < B; i++)
			tour[pos++] = tmp[i];
	}

	public static class IntegerTriple implements Comparable<IntegerTriple> {
		int _first;
		int _second;
		int _third;

		public IntegerTriple(int f, int s, int t) {
			_first = f;
			_second = s;
			_third = t;
		}

		public int compareTo(IntegerTriple o) {
			if (first() != o.first())
				return this.first() - o.first();
			else if (second() != o.second())
				return this.second() - o.second();
			else
				return this.third() - o.third();
		}

		int first() {
			return _first;
		}

		int second() {
			return _second;
		}

		int third() {
			return _third;
		}
	}

	public static class MyScanner {
		private static BufferedReader br;
		private String line;
		private StringTokenizer st;

		public MyScanner() {
			br = new BufferedReader(new InputStreamReader(System.in));
		}

		String nextLine() {
			try {
				line = br.readLine();
			} catch (IOException e) {
				e.printStackTrace();
			}

			return line;
		}

		String next() {
			while (st == null || !st.hasMoreElements()) {
				nextLine();

				if (line == null) {
					return null;
				}

				st = new StringTokenizer(line);
			}

			return st.nextToken();
		}

		int nextInt() {
			try {
				return Integer.parseInt(next());
			} catch (Exception e) {
				return -1;
			}
		}

		long nextLong() {
			try {
				return Long.parseLong(next());
			} catch (Exception e) {
				return -1;
			}
		}

		double nextDouble() {
			try {
				return Double.parseDouble(next());
			} catch (Exception e) {
				return -1;
			}
		}
	}
}