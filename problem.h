/*
 * problem.h
 *
 *  Created on: 2016/02/08
 *      Author: aman
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

namespace allocation_optimize_NS {

class Problem {
public:
	class DC {
	private:
		static int obj_cnt;
		explicit DC(int cmp) { // 不要なコンストラクタ
			idx = obj_cnt;
			cmp_c = cmp;
			capa = 0;
			obj_cnt++;
		}
	public:
		int idx; // インデックスはオブジェクトの生成順に番号が振られる
		int cmp_c; // computation cost
		int capa;  // capacity
		std::vector<std::pair<int, DC*>> adjacentNodes;

		DC(int cmp, int c) {
			idx = obj_cnt;
			cmp_c = cmp;
			capa = c;
			obj_cnt++;
		}
		~DC(){
			obj_cnt--;
		}
		static int cnt() { return obj_cnt; }
		bool operator<(const DC &right) const { return this->idx < right.idx; }
	};

	// Work flow task
	class Task {
	private:
		static int obj_cnt;
		/*
		Task() { // 不要なコンストラクタ
			idx = obj_cnt;
			cmp_r = 0;
			cmm_r = 0;
			assignNode = nullptr;
			obj_cnt++;
		}
		 */
	public:
		int idx, cmp_r, cmm_r; // computation_requiremt & communication_requirement
		DC* assignNode;
		Task(int cmp, int cmm) {
			idx = obj_cnt;
			cmp_r = cmp;
			cmm_r = cmm;
			assignNode = nullptr;
			obj_cnt++;
		}
		~Task(){
			assignNode = nullptr;
			obj_cnt--;
		}
		// Task(const Task &obj) { }
		static int cnt() { return obj_cnt; }

		bool operator<(const Task &right) const { return this->idx < right.idx; }
	};

	// 得られた解のコストの値と配置した内訳を表示
	class Solution {
	private:
		static int obj_cnt;
	protected:
		int resultCommCost;
		int resultCompCost;
		std::list < std::pair < Task*, DC* >> allocationList; // 各タスクの配置した内訳
	public:
		Solution() : resultCommCost{ 0 }, resultCompCost{ 0 } { }
		Solution(int com, int cmp) : resultCommCost{ com }, resultCompCost{ cmp } { }
		~Solution() {
			obj_cnt--;
		}
		/*
			std::pair<int, int> resultCost() {
			return{ resultCommCost, resultCompCost };
			}
		 */
		void manipCommCost(int comm) { resultCommCost += comm; }
		void manipCompCost(int comp) { resultCompCost += comp; }
		int cntAllocatedNodes() { return allocationList.size(); }
		void regAllocation(Task *task, DC *dc) {
			std::cout << "Task" << task->idx << " = {" << task->cmp_r << ", " << task->cmm_r << "}" << std::endl;
		//	std::cout << task->idx << " vs " << dc->idx << std::endl;
			//allocationList.push_back(std::map<Task*,DC*>::value_type(t, d));
			if(dc->capa >= task->cmp_r) {
				std::cout << "DC" << dc->idx << " = {" << dc->capa << ", " << dc->cmp_c << "}";
				allocationList.push_back({ task, dc });
				dc->capa -= task->cmp_r;
				std::cout << " => " << "{" << dc->capa << ", " << dc->cmp_c << "}"  << std::endl;
			}
		}
		int getResultCost() { return resultCommCost + resultCompCost; }

		//void printAllocations(std::vector<Task*> ti) {
		void print() {
			printResultCost();
			printAllocation();
		}
		void printResultCost() {
			std::cout << "Total computation cost = " << resultCompCost << std::endl;
			std::cout << "Total cost = " << getResultCost() << std::endl;
		}
		void printAllocation() {
			std::cout << "Size of allocation list: " << allocationList.size() << std::endl;
			//auto std::list<std::pair<Task*, DC*>>::iterator it;
			//auto it = allocationList.begin();
			for (auto pair : allocationList) {
				//for (uint32_t i = 0; i < allocationList.size(); ++i)
				//std::cout << "Task" << allocationList << ": " << "DC" << allocationList[i]-> << std::endl;
				std::cout << "Task" << pair.first->idx << " : " << "DC" << pair.second->idx << std::endl;
			}
		}
	};

	// 全てのタスクがDCに割り当てられたら(=タスク配置がすべて終了したら)
	//bool allAssigned(const std::vector<Task*>& ti) {
	bool allAssigned() {
		bool flag = true;
		for (auto task : ti)
			if (task->assignNode == nullptr)
				flag = false;
		return flag;
	}
private:
	int N{}, M{}; // Task amount N & DC node amount M
	Problem::DC *maxDegree;

	// Set of DC & Task
public:
	std::vector<Task*> ti; // workflow
	std::vector<Problem::DC*> dj; // set of data centers (=network topology)
	std::vector<Solution*> solutions; // set of data centers (=network topology)
	// N,Mを指定しなければプリセットを用いて初期化する
	// If you had not set the values of N and M, to initialize by preset status.

	Problem() {
		// Task list initialize => Cmp_R, Cmm_R
		ti = { new Task(4, 10), new Task(7, 15), new Task(8, 5), new Task(5, 10) };
		// DC nodes initialize => Cmp_C, Capacity
		dj = { new DC(10, 5), new DC(14, 7), new DC(8, 3), new DC(15, 8), new DC(20, 10), new DC(9, 8) };

		N = ti.size(); M = dj.size();

		// Adjacent relations definition
		dj[0]->adjacentNodes = {
				{ 3, dj[1] },
				{ 2, dj[3] },
				{ 4, dj[4] }
		};
		dj[1]->adjacentNodes = {
				{ 3, dj[0] },
				{ 3, dj[2] },
				{ 4, dj[3] }
		};
		dj[2]->adjacentNodes = {
				{ 3, dj[1] },
				{ 6, dj[3] },
				{ 3, dj[4] },
				{ 1, dj[5] }
		};
		dj[3]->adjacentNodes = {
				{ 2, dj[0] },
				{ 4, dj[1] },
				{ 6, dj[2] }
		};
		dj[4]->adjacentNodes = {
				{ 5, dj[0] },
				{ 3, dj[2] },
				{ 2, dj[5] }
		};
		dj[5]->adjacentNodes = {
				{ 1, dj[2] },
				{ 2, dj[4] }
		};
	}

	// int型変数indexで指定された番号のDCノードを求める
	// 注意：djの初期化が完全に終了してからでないとsegmentation faultを引き起こす可能性がある
	DC* getDCIdx(int idx) {
		return dj[idx];
	}

	// 最も多くのノードとつながる次数のノードを求める
	DC* getMaxDegreeNode() {
		return [=](std::vector<DC*> adjNodes){
			unsigned int max = std::numeric_limits<unsigned int>::min();
			DC* tmp;
			for (auto dst : adjNodes) {
				//std::cout << d->adjacentNodes.size() << std::endl;
				if (dst->adjacentNodes.size() >= max) {
					max = dst->adjacentNodes.size();
					tmp = dst;
				}
			}
			return tmp;
		}(dj);
	}

	// 読み込んだjsonファイルからProblemクラスを設定するコンストラクタ
	Problem(std::string file) {
		std::cout << "Input filename: " << file << std::endl << std::endl;
		std::ifstream stream(file);
		if(!stream.is_open()) exit(-1);
		picojson::value root; // json file root
		stream >> root;

		// 入力するファイルはプログラムと同じ階層にある
		picojson::object problem = root.get<picojson::object>()["Problem"].get<picojson::object>();
		picojson::array task = problem["Task"].get<picojson::array>();
		for(uint32_t i = 0; i < task.size(); ++i) {
			picojson::object obj = task[i].get<picojson::object>();
			int cmp_r = static_cast<int>(obj["cmp_r"].get<double>());
			int cmm_r = static_cast<int>(obj["cmm_r"].get<double>());
			ti.push_back(new Task(cmp_r, cmm_r));
		}

		picojson::array dc = problem["DC"].get<picojson::array>();
		for(uint32_t j = 0; j < dc.size(); ++j) {
			picojson::object obj = dc[j].get<picojson::object>();
			int cmp_c = static_cast<int>(obj["cmp_c"].get<double>());
			int capa = static_cast<int>(obj["capa"].get<double>());
			dj.push_back(new DC(cmp_c, capa));
		}

		// いったんDCのコンテナであるdjにノード情報を格納してから出ないと
		// 互いの隣接情報を設定できないので分けてある
		for(uint32_t j = 0; j < dc.size(); ++j) {
			picojson::array adjacentNodes = dc[j].get<picojson::object>()["adjacentNodes"].get<picojson::array>();
			for(uint32_t idx = 0; idx < adjacentNodes.size(); ++idx){
				picojson::object obj = adjacentNodes[idx].get<picojson::object>();
				int cmm_c = static_cast<int>(obj["cmm_c"].get<double>());
				int adjIdx = static_cast<int>(obj["adjIdx"].get<double>());
				std::pair<int, Problem::DC*> tmp = {cmm_c, getDCIdx(adjIdx) };
				dj[j]->adjacentNodes.push_back({ cmm_c, getDCIdx(adjIdx) });
			}
		}
		N = ti.size(); M = dj.size();
	}

	~Problem() {
		for(auto task : ti) {
			delete task;
		}
		for(auto dc: dj) {
			delete dc;
		}
		for(auto s: solutions) {
			delete s;
		}
	}



	// Task数NとDCノード数Mが与えられた場合、ランダムに生成したcmp_r/c, cmm_r/cの値を使って初期化する
	// Input the amout of Task N and DC node M,
	// then each (Cmp/Cmm)_(requirement/cost) values will initialized by random integer.
	/*
		void init(int n, int m) {
		// this lines used for random number library
		std::random_device rnd;
		std::mt19937 mt(rnd());
		std::uniform_int_distribution<> range(1,100); // random number based on distribution
		N = n; M = m;

		for(int i = 0; i < N; ++i) ti.push_back(new Task(range(mt),range(mt)));
		for(int j = 0; j < M; ++j) dj.push_back(new DC(range(mt)/10));

		edge.resize(N);
		for(int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
		edge[i].push_back(range(mt)/10); // communication_cost
		}
		}
		path.resize(N);
		for(int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
		path[i].push_back(
		std::tuple<int,float,int>(
		dj[j]->capa,
		static_cast<float>(dj[j]->cmp_c+edge[i][j]),
		0
	}
	 */

	// Task数NとDCノード数Mが与えられた場合、ランダムに生成したcmp_r/c, cmm_r/cの値を使って初期化する
	// Input the amout of Task N and DC node M,
	// then each (Cmp/Cmm)_(requirement/cost) values will initialized by random integer.

	/*
		void init(int n, int m) {
		// this lines used for random number library
		std::random_device rnd;
		std::mt19937 mt(rnd());
		std::uniform_int_distribution<> range(1,100); // random number based on distribution
		N = n; M = m;

		for(int i = 0; i < N; ++i) ti.push_back(new Task(range(mt),range(mt)));
		for(int j = 0; j < M; ++j) dj.push_back(new DC(range(mt)/10));

		edge.resize(N);
		for(int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
		edge[i].push_back(range(mt)/10); // communication_cost
		}
		}
		path.resize(N);
		for(int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
		path[i].push_back(
		std::tuple<int,float,int>(
		dj[j]->capa,
		static_cast<float>(dj[j]->cmp_c+edge[i][j]),
		0
		)
		);
		}
		}
		}
	 */

	void print() {
		std::cout << "Task amount: " << Problem::Task::cnt()
		<< ", DC nodes amount: " << Problem::DC::cnt() << std::endl;

		std::cout << "computation_req, communication_req" << std::endl;
		for (int i = 0; i < N; ++i)
			std::cout << "task" << i << ":" << ti[i]->cmp_r << ", " << ti[i]->cmm_r << std::endl;

		std::cout << std::endl << "computation_costs, capacity, adjacency info" << std::endl;
		for (auto dc : dj) {
			std::cout << "DC" << dc->idx << ":" << dc->cmp_c << ", " << dc->capa << " -> { ";
			for (uint32_t j = 0; j < dc->adjacentNodes.size(); ++j)
				//std::cout << "DC" << dc->adjacentNodes[j].second->idx << "(" << dc->adjacentNodes[j].first << ") ";
				std::cout << "DC" << dc->adjacentNodes[j].second->idx << "(" << dc->adjacentNodes[j].first << ") ";
			std::cout << "}" << std::endl;
		}
		std::cout << std::endl;
		std::cout << "max degree node: DC" << Problem::getMaxDegreeNode()->idx << " has " << Problem::getMaxDegreeNode()->adjacentNodes.size() << " nodes." << std::endl;
	}

	// 要求されたcmp_rに対してキャパシティーの余裕があり且つ最小のコストでたどり着けるノードを探索する
	DC* minCmpCost(int cmp_r) {//, std::vector<DC*> dj) {
		auto min = std::numeric_limits<int>::max();
		int minIdx{};

		for (auto dc : dj) {
			if (dc->cmp_c <= min && dc->capa >= cmp_r) {
				min = dc->cmp_c;
				minIdx = dc->idx;
			}
			// ↓間違い、すべてのタスクを配置できることが前提条件で与えられている
			// else // 配置できるノードが一つもなかったら
		}
		std::cout << "cmp_r= " << cmp_r << std::endl;
		std::cout << "minIdx DC[" << minIdx << "], cmp_c= " << dj[minIdx]->cmp_c << ", capa= " << dj[minIdx]->capa << std::endl;
		return dj[minIdx];
	}

	class Timer {
	private:
		std::chrono::time_point<std::chrono::system_clock, std::chrono::system_clock::duration> start, end;
	public:
		void timerStart() {
			start = std::chrono::system_clock::now();
		}
		void timerEnd(bool flag) {
			end = std::chrono::system_clock::now();
			if (flag) {
				auto elapse = end - start;
				std::cout << "duration time : " << std::chrono::duration_cast<std::chrono::milliseconds>(elapse).count() << " msec" << std::endl;
			}
		}

	};
	//Problem::Solution*
	//void optimalAllocationPattern(std::list<Problem::Solution> solutions) {
	void optimalAllocationPattern() {
		int minCost = std::numeric_limits<int>::max();
		Problem::Solution *min_ptr = nullptr;

		for (auto s : this->solutions) {
			if (s->getResultCost() <= minCost) {
				std::cout << s->getResultCost() << std::endl;
				minCost = s->getResultCost();
				min_ptr = s;
			}
		}
		if (min_ptr != nullptr) {
			std::cout << "optimal cost value: " << min_ptr->getResultCost() << std::endl;
			//printAllocations();
		}
	}

}; // end of class Problem difinition
int Problem::DC::obj_cnt;
int Problem::Task::obj_cnt;
int Problem::Solution::obj_cnt;

// 総当たりによる最適解の探索
// 探索範囲によっては計算が終了しないが絶対に最適解を見つけられる
// Searching of optimal solution by brute force method
/*
std::list<Problem::Solution> brute_force(Problem *p) {
	Problem initialState = *p;
	Problem tmpProblem = *p;
	std::list<Problem::Solution> solutions; // 解の集合

	uint32_t cnt{ 0 };
	//for (cnt = 0; cnt < pow(p->dj.size(), p->ti.size()); ++cnt) { // 全パターンはmΠn(m=dc amount,n=task amount)+alpha

	//while(!tmpProblem.allAssigned()) {
		Problem::Solution tmpSolution;
		//int cmpC{ 0 }, cmmC{ 0 }, minCost = std::numeric_limits<int>::max();
		//		p->init(); // グラフ状況の初期化
		for (auto task : tmpProblem.ti) {
			int assignStatus = -1;
			for(auto dc : tmpProblem.dj) {
			//if (p->allAssigned()) break;
			//Problem::DC *tmpMinimumNode = p->minCmpCost(task->cmp_r);
				assignStatus = tmpSolution.regAllocation(task, dc);
				if(assignStatus==0) {
					std::cout << "Task" << task->idx << " assigned to DC" << dc->idx << std::endl;
					break;
				}
			}
			if(assignStatus==1) break; // もう配置できない
		}

		//noneAllocatableNodes()
	//	if([=]()->bool{}()) break;

		//solutions.push_back(tmpSolution);
		//tmpProblem = *p;
	//}
	std::cout << "Patterns of solution: " << cnt << std::endl;
	//tmpProblem.optimalAllocationPattern();
	//solutions.push_back(tmpSolution);
	p->optimalAllocationPattern();

	return solutions;
}
*/

// 深さ優先探索で与えられた2つのノード間の経路を調べる
// この探索で得られた経路は最小cmm_cとは限らないので注意
std::list<Problem::DC *> dfsSearch(Problem *p, Problem::DC *src, Problem::DC *dst) {
	std::vector<std::list<Problem::DC *>> allPaths{};
	std::list<Problem::DC *> path{};
	std::stack<Problem::DC *> q;

	q.push(src); // q0 push
	// std::cout << "src: " << src->idx << ", dst: " << dst->idx << std::endl;

	std::list<Problem::DC *> visitedMemory;
	auto isVisited = [=](Problem::DC *v, std::list<Problem::DC *> mem) {
		for(auto visited : mem) {
			if(visited->idx == v->idx) return true;
		}
		return false;
	};

	while(!q.empty()) {
		Problem::DC *v = q.top();
		path.push_back(v);
		visitedMemory.push_back(v);
		q.pop();

		if(v->idx == dst->idx) {
			//allPaths.push_back(visitedMemory);
			//q.pop();
			break;
		}
		else {
			for(auto adj : v->adjacentNodes) {
				// 未訪問ノードのうち、最もコストの小さいもの
				//if(!isVisited(adj.second, visitedMemory)) q.push(adj.second);
				// この実装では各隣接ノードを調べる段階で最もコストの低いものしか選べないので
				// greedyな探索である
				std::list<std::pair<int, Problem::DC *>> nonVisitedAdj{};
				for(auto dc : v->adjacentNodes) {
				if(!isVisited(dc.second, visitedMemory)) nonVisitedAdj.push_back(dc);
				}

				Problem::DC *minCmmAdj  = [=](){
					int min = std::numeric_limits<int>::max();
					Problem::DC *minAdj = nullptr;
					for(auto adjInfo : nonVisitedAdj) {
						if(min >= adjInfo.first) {
							min = adjInfo.first;
							minAdj = adjInfo.second;
						}
					}
					//std::cout << minAdj->idx << std::endl;
					return minAdj;
				}();

				//q.push(adj.second);
				q.push(minCmmAdj);
			}
		}
		v = nullptr;
	}

	for(auto dc : path) {
		std::cout << dc->idx;
		if(dc != dst) std::cout << "->";
		else std::cout << std::endl;
	}
	std::cout << std::endl;

	return path;
}

std::list<Problem::Solution> brute_force(Problem *p) {
	std::list<Problem::Solution> solutions; // 解の集合

	std::list<Problem::DC *> path{};
	//path = dfsSearch(p, p->dj[2], p->dj[0]);

	/*
	uint32_t cnt{ 0 };
	for (cnt = 0; cnt < pow(p->dj.size(), p->ti.size()); ++cnt) { // 全パターンは最大でmΠn(m=dc amount,n=task amount)+alpha分存在する
		int cmpC{ 0 }, cmmC{ 0 }, minCost = std::numeric_limits<int>::max();
		//		p->init(); // グラフ状況の初期化
		Problem::Solution tmpSolution;
		for (auto task : p->ti) {
			if (p->allAssigned()) break;
			Problem::DC *tmpMinimumNode = p->minCmpCost(task->cmp_r);
			std::cout << "Temporary: Search() End" << std::endl;
			tmpSolution.regAllocation(task, tmpMinimumNode);
		}
		solutions.push_back(tmpSolution);
	}
	std::cout << "Patterns of solution: " << cnt << std::endl;
	p->optimalAllocationPattern();
	*/
	while(!p->allAssigned()) {

	}

	return solutions;
}


// greedy method
// 貪欲法によるcmp_cによるタスク配置
// 未完成：その時点で選べる最小コストのノードを選ぶようになっていない
Problem::Solution* greedy(Problem *p) {
	Problem::Solution* solution = new Problem::Solution();
	int cmpC{ 0 }, cmmC{ 0 };
	int calc_steps = 0;

		for(auto task : p->ti) {
			// task[i]のcmp_costを許容できるcapaを持つDCのうち、cmpコストが最小のもの
			Problem::DC *tmpMinCostNode = [=]{
				int minCost = std::numeric_limits<int>::max();
				Problem::DC *minCmpCostNode{};
				for(auto dc : p->dj) {
					if(dc->capa >= task->cmp_r && dc->cmp_c <= minCost) {
						minCost = dc->cmp_c;
						minCmpCostNode = dc;
					}
				}
				return minCmpCostNode;
			}();
			solution->regAllocation(task, tmpMinCostNode);
			if(p->allAssigned()) break;
		}

	solution->printAllocation();
	std::cout << std::endl;

	return solution;
}


/*

// greedy method
std::pair<int, int> greedy(Problem *p) {
	int cmpC{ 0 }, cmmC{ 0 }, minCost = std::numeric_limits<int>::max();
	int calc_steps = 0;
	Problem::DC* currPosition = nullptr; // 最後に配置したノードの番号

	// 未完成：その時点で選べる最小コストのノードを選ぶようになっていない

	// 最初に配置するスタート地点の設定
	currPosition = p->minCmpCost(0);
	// std::cout << currPosition->idx << std::endl;

	// 隣接するなかで最小のノードを選ぶ

	std::cout << "Allocation result: " << std::endl;
	for (auto task : p->ti) {
		if (p->allAssigned()) break;
		for (auto dc : p->dj) {
			// 配置できるのなら
			if ((dc->capa >= task->cmp_r) && (task->assignNode == nullptr)) {

				cmpC += dc->cmp_c*task->cmp_r;
				dc->capa -= task->cmp_r;
				task->assignNode = dc;
				currPosition = dc;

				std::cout << "task" << task->idx << " -> " << "DC" << dc->idx << std::endl;
			}
		}
	}
	return std::pair<int, int>(cmpC + cmmC, calc_steps);
}
*/

} /* allocation_optimize_NS */

#endif /* PROBLEM_H_ */
