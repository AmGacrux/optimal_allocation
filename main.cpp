#include <cstdlib>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <list>
#include <tuple>
#include <random>
#include <vector>
#include <stack>

// Allocation problem and algorithms (under construction)
// 
// Compile Requirement: g++ -std=c++11 allocation.cpp
// 

namespace allocation_optimize_NS {

	// Data center node
	class DC {
	private:
		static int obj_cnt;
	public:
		int idx; // computation cost
		int cmp_c; // computation cost
		int capa;  // capacity
		std::vector<std::pair<int, DC*>> adjacentNodes;

		explicit DC(int cmp) {
			idx = obj_cnt;
			cmp_c = cmp;
			capa = 0;
			obj_cnt++;
		}
		DC(int cmp, int c) {
			idx = obj_cnt;
			cmp_c = cmp;
			capa = c;
			obj_cnt++;
		}
		static int cnt() { return obj_cnt; }
		bool operator<(const DC &right) const { return this->idx < right.idx; }
	};
	int DC::obj_cnt{ 0 };

	/*
	// Data center node
	class Path {
	private:
	static int obj_cnt;
	public:
	int cmm_c; // computation cost
	std::pair<DC*,DC*> endPoints;
	explicit Path(int cmm) {
	cmm_c = cmm;
	obj_cnt++;
	}
	Path(int cmm, const DC& d1, const DC& d2) {
	cmm_c = cmm;
	endPoints.first = const_cast<DC*>(&d1);
	endPoints.second = const_cast<DC*>(&d2);
	obj_cnt++;
	}
	static int cnt() { return obj_cnt; }
	};
	int Path::obj_cnt;
	*/

	// Work flow task
	class Task {
	private:
		static int obj_cnt;
	public:
		int idx, cmp_r, cmm_r; // computation_requiremt & communication_requirement
		DC* assignNode;
		Task() {
			idx = obj_cnt;
			cmp_r = 0;
			cmm_r = 0;
			obj_cnt++;
		}
		Task(int cmp, int cmm) {
			idx = obj_cnt;
			cmp_r = cmp;
			cmm_r = cmm;
			assignNode = nullptr;
			obj_cnt++;
		}
		Task(const Task &obj) { }
		static int cnt() { return obj_cnt; }
		bool operator<(const Task &right) const { return this->idx < right.idx; }
	};
	int Task::obj_cnt{ 0 };

	// 全てのタスクがDCに割り当てられたら(=タスク配置がすべて終了したら)
	bool allAssigned(const std::vector<Task*>& ti) {
		bool flag = true;
		for (auto task : ti) if (task->assignNode == nullptr) flag = false;
		return flag;
	}

	// Set of DC & Task
	int N, M; // Task amount N & DC node amount M
	int maxDegree{};
	std::vector<Task*> ti; // workflow
	std::vector<DC*> dj; // set of data centers (=network topology)

	//	std::vector<std::vector<float>> edge;
	//	std::vector<std::vector<std::tuple<int,float,int>>> path; // (capacity, cost=cmp_c+cmm_c, flow)
	// 	std::vector<std::pair<int,int>> alloc;

	// N,Mを指定しなければプリセットを用いて初期化する
	// If you had not set the values of N and M, to initialize by preset status.
	void init() {
		// Task list initialize
		ti.push_back(new Task(4, 10)); // Cmp_R, Cmm_R
		ti.push_back(new Task(7, 15));
		ti.push_back(new Task(8, 5));
		ti.push_back(new Task(5, 10));

		// DC nodes initialize
		dj.push_back(new DC(10, 5)); // Cmp_C, Capacity
		dj.push_back(new DC(14, 7));
		dj.push_back(new DC(8, 3));
		dj.push_back(new DC(15, 8));
		dj.push_back(new DC(20, 10));
		dj.push_back(new DC(9, 8));

		N = ti.size(); M = dj.size();

		// Adjacent relations definition
		dj[0]->adjacentNodes = {
			{ 3, dj[1] },
			{ 2, dj[3] },
			{ 4, dj[4] }
		};
		/*
		// The above lines same these
		dj[0]->adjacentNodes.push_back({3, dj[1]}); // initializer list (it behave a std::pair<int,DC*>)
		dj[0]->adjacentNodes.push_back({2, dj[3]});
		dj[0]->adjacentNodes.push_back({5, dj[4]});
		*/

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
		maxDegree = [=](std::vector<DC*> adj){
			unsigned int max = std::numeric_limits<unsigned int>::min();
			for (auto d : adj) {
				//std::cout << d->adjacentNodes.size() << std::endl;
				if (d->adjacentNodes.size() >= max) max = d->adjacentNodes.size();
			}
			return max;
		}(dj);
		//std::cout << "DJ:" << maxDegree << std::endl;

		/*
		edge.resize(N);
		edge[0].push_back(0.1f); // Cmm_c
		edge[0].push_back(0.2f);
		edge[0].push_back(0.2f);
		edge[0].push_back(0.8f);
		edge[1].push_back(0.1f);
		edge[1].push_back(0.3f);
		edge[1].push_back(0.3f);
		edge[1].push_back(0.3f);
		edge[2].push_back(0.2f);
		edge[2].push_back(0.2f);
		edge[2].push_back(0.3f);
		edge[2].push_back(0.5f);
		*/

		/*
		path.resize(N);
		for(int i = 0; i < N; ++i) {
		for(int j = 0; j < M; ++j) {
		path[i].push_back(
		std::tuple<int, float, int>(
		dj[j]->capa,
		static_cast<float>(dj[j]->cmp_c+edge[i][j]),
		0
		)
		);
		}
		}
		*/
	}

	// 得られた解のコストの値と配置した内訳を表示
	class Solution {
	private:
		static int obj_cnt;
		int resultCommCost; int resultCompCost;
		std::list < std::pair < Task*, DC* >> allocationList; // 各タスクの配置した内訳
	public:
		Solution() : resultCommCost{ 0 }, resultCompCost{ 0 } { }
		Solution(int com, int cmp) : resultCommCost{ com }, resultCompCost{ cmp } { }
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
		void printAllocations(std::vector<Task*> ti, std::vector<DC*> dj) {
			for (auto task : ti) std::cout << task->assignNode << std::endl;
		}
	};
	int Solution::obj_cnt{ 0 };

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
	}

	// 		for(auto &a : dj) std::cout << a->cmp_c << " ";
	// 		std::cout << std::endl;
	/*
	std::cout << "communication_costs" << std::endl;
	for(int i = 0; i < N; ++i) {
	for(int j = 0; j < M; ++j) {
	std::cout << "t" << i << "->d" << j << "=";
	std::cout << edge[i][j] << " ";
	}
	std::cout << std::endl;
	}

	std::cout << std::endl << "path(capacity, integrated cost, flow)" << std::endl;
	for(int i = 0; i < N; ++i) {
	for(int j = 0; j < M; ++j) {
	std::cout << "(" << std::get<0>(path[i][j]) << ", "
	<< std::get<1>(path[i][j])<< ", "
	<< std::get<2>(path[i][j]) << ")" << " ";
	}
	std::cout << std::endl;
	}
	*/

	/*
	bool allocate_task(DC* d, Task* t) {
	if (d->capa > t->cmp_r) {
	// 			d->assigned_tasks.push_back(t);
	t->assignNode = d;
	d->capa -= t->cmp_r;
	return true;
	}
	else return false;
	}
	*/

	// 要求されたcmp_rに対してキャパシティーの余裕があり且つ最小のコストでたどり着けるノードを探索する
	DC* minCmpCost(int cmp_r, std::vector<DC*>& dj) {
		auto min = std::numeric_limits<int>::max();
		int minIdx{};

		for (auto dc : dj) {
			if (dc->cmp_c <= min && dc->capa >= cmp_r) {
				min = dc->cmp_c;
				minIdx = dc->idx;
			}
		}
		return dj[minIdx];
	}

	// 総当たりによる最適解の探索
	// 探索範囲によっては計算が終了しないが絶対に最適解を見つけられる
	// Searching of optimal solution by brute force method
	std::list<Solution>& brute_force() {
		int cmpC{ 0 }, cmmC{ 0 }, minCost = std::numeric_limits<int>::max();
		//int idx, calc_steps = 0;
		std::vector<Solution> solutions;

		Solution tmpSolution;

		return solutions;
	}

	// greedy method
	std::pair<float, int> greedy() {
		int cmpC{ 0 }, cmmC{ 0 }, minCost = std::numeric_limits<int>::max();
		int tmp, idx, calc_steps = 0;
		DC* currPosition = nullptr; // 最後に配置したノードの番号

		// 未完成：その時点で選べる最小コストのノードを選ぶようになっていない

		// 最初に配置するスタート地点の設定
		currPosition = minCmpCost(0, dj);
		// std::cout << currPosition->idx << std::endl;

		std::cout << "Allocation result: " << std::endl;
		for (auto task : ti) {
			if (allAssigned(ti)) break;
			for (auto dc : dj) {
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
			std::cout << "duration time : "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(elapse).count()
				<< " msec" << std::endl;
		}

	}
};

int main(int argc, char** argv) {
	namespace OptNS = allocation_optimize_NS;

	//if(argc > 2) allocation_optimize_NS::init(atoi(argv[1]), atoi(argv[2]));
	//else if(argc == 1) 
	allocation_optimize_NS::init();

	allocation_optimize_NS::print();
	std::cout << std::endl << "Task amount: " << OptNS::Task::cnt() << ", DC nodes amount: " << OptNS::DC::cnt() << std::endl;

	Timer timer;

	std::cout << "================================" << std::endl;
	/*
	std::cout << "Using greedy method:" << std::endl;
	timer.timerStart();
	std::pair<int, int> ans = allocation_optimize_NS::greedy();
	timer.timerEnd(true);
	//*/

	std::cout << "Using brute force method:" << std::endl;
	timer.timerStart();
	std::list<allocation_optimize_NS::Solution> solutions = allocation_optimize_NS::brute_force();
	timer.timerEnd(true);

	//std::cout << "result: " << ans.first << std::endl;
	//std::cout << "steps: " << ans.second << std::endl;
	int N{};
	std::cin >> N;

	return int32_t(0);
}