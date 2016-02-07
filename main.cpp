#include <chrono>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <stack>
#include <string>
#include <tuple>
#include <random>
#include <vector>

#include <picojson.h> // for read a json file
#include <fstream>

// Allocation problem and algorithms (under construction)
// 
// Compile Requirement: g++ -std=c++11 allocation.cpp
// 

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
		void regAllocation(Task *t_ptr, DC *d_ptr) {
			//allocationList.push_back(std::map<Task*,DC*>::value_type(t, d));
			allocationList.push_back({ t_ptr, d_ptr });
		}
		int getResultCost() { return resultCommCost + resultCompCost; }

		//void printAllocations(std::vector<Task*> ti) {
		void print() {
			printResultCost();
			printAllocations();
		}
		void printResultCost() {
			std::cout << "Total computation cost = " << resultCompCost << std::endl;
			std::cout << "Total cost = " << getResultCost() << std::endl;
		}
		void printAllocations() {
			//auto std::list<std::pair<Task*, DC*>>::iterator it;
			//auto it = allocationList.begin();
			for (auto idx : allocationList) {
				//for (uint32_t i = 0; i < allocationList.size(); ++i)
				//std::cout << "Task" << allocationList << ": " << "DC" << allocationList[i]-> << std::endl;
				std::cout << "Task" << idx.first << ": " << "DC" << idx.second << std::endl;
			}
		}
	};

	// 全てのタスクがDCに割り当てられたら(=タスク配置がすべて終了したら)
	//bool allAssigned(const std::vector<Task*>& ti) {
	bool allAssigned() {
		bool flag = true;
		for (auto task : ti) if (task->assignNode == nullptr) flag = false;
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

		// 最も多くのノードとつながる次数を求める
		maxDegree = [=](std::vector<DC*> adjNodes){
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
		std::cout << "max degree node: DC" << maxDegree->idx << " has " << maxDegree->adjacentNodes.size() << " nodes." << std::endl;
	}

	// 未完成
	// Todo: 読み込んだjsonファイルからProblemクラスを設定するコンストラクタ
	Problem(std::string file) {
		std::cout << "Input filename: " << file << std::endl;
		picojson::value root;
		// 入力するファイルはプログラムと同じ階層にある
		std::ifstream stream(file);
		stream >> root;
		picojson::object problem = root.get<picojson::object>()["Problem"].get<picojson::object>();
		//std::cout << problem["Graph"].get<std::string>() << std::endl;
		//picojson::array ti = task["idx"].get<picojson::array>();
		picojson::array ti = problem["Task"].get<picojson::array>();
		  // arrayはstd::vector<picojson::value>なのでrange-based forが使える!
		  for (picojson::value task: ti) {
		    std::cout << task.get("idx") << std::endl;
		  }
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
		std::cout << std::endl << "Task amount: " << Problem::Task::cnt()
		<< ", DC nodes amount: " << Problem::DC::cnt() << std::endl;

		std::cout << "computation_req, communication_req" << std::endl;
		for (int i = 0; i < N; ++i)
			std::cout << "task" << i << ":" << ti[i]->cmp_r << ", " << ti[i]->cmm_r << std::endl;

		std::cout << std::endl << "computation_costs, capacity, adjacency info" << std::endl;
		for (auto dc : dj) {
			std::cout << "DC" << dc->idx << ":" << dc->cmp_c << ", " << dc->capa << " -> { ";
			for (uint32_t j = 0; j < dc->adjacentNodes.size(); ++j)
				std::cout << "DC" << dc->adjacentNodes[j].second->idx << " ";
			std::cout << "}" << std::endl;
		}
		std::cout << std::endl;
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
std::list<Problem::Solution> brute_force(Problem *p) {
	std::list<Problem::Solution> solutions; // 解の集合

	uint32_t cnt{ 0 };
	for (cnt = 0; cnt < pow(p->dj.size(), p->ti.size()); ++cnt) { // 全パターンはmΠn(m=dc amount,n=task amount)
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

	return solutions;
}

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

} // end of namespace allocation_optimize_NS

int main(int argc, char** argv) {
	//if(argc > 2) allocation_optimize_NS::init(atoi(argv[1]), atoi(argv[2]));
	//else if(argc == 1) 

	//if(argc == 2) std::string file{argv[1]};
	std::string file{"dataset_2.json"};

	//std::istringstream cmd("");
	std::string cmd;
	while (true) {
allocation_optimize_NS::Problem *p{nullptr};

		if(!file.empty()) p = new allocation_optimize_NS::Problem(file);
		else if(argc < 1) p = new allocation_optimize_NS::Problem();
		p->print();

		std::cout << "Choose a search algorithm. (g:Greedly, b:Brute force, ... , q:Quit this program) > ";
		std::cin >> cmd;

		if (cmd == "q") {
			std::cout << "end." << std::endl;
			delete p;
			return uint32_t(0);
		}

		//std::cout << "================================" << std::endl;
		///*
		else if(cmd == "g") {
			allocation_optimize_NS::Problem::Timer timer;
			std::cout << "Using greedy method:" << std::endl;
			timer.timerStart();
			std::pair<int, int> ans = allocation_optimize_NS::greedy(p);
			timer.timerEnd(true);

			std::cout << "greedy optimization is end." << std::endl;
			delete p;

			continue;

			//std::cin >> cmd;
			//return uint32_t(0);
		}
		else if(cmd == "b") {
			allocation_optimize_NS::Problem::Timer timer;
			std::cout << "Using brute force method:" << std::endl;
			//std::list<allocation_optimize_NS::Problem::Solution> solutions = allocation_optimize_NS::brute_force(p);
			timer.timerEnd(true);

			delete p;
			continue;

			//std::cin >> cmd;
			//return uint32_t(0);
		}
		else {
			std::cout << "Invalid command input..." << std::endl;

			delete p;
			continue;
		}
		//std::cout << "result: " << ans.first << std::endl;
		//std::cout << "steps: " << ans.second << std::endl;
	}

	return int32_t(0);
}
