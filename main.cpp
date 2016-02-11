#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <stack>
#include <string>
#include <tuple>
#include <random>
#include <vector>
#include <queue>

#include "picojson.h" // for read a json file

#include "problem.h" // 自作クラスのヘッダ

// Allocation problem and algorithms (under construction)
//
// Compile Requirement: g++ -std=c++11 allocation.cpp
//

int main(int argc, char** argv) {
	//if(argc > 2) allocation_optimize_NS::init(atoi(argv[1]), atoi(argv[2]));
	//else if(argc == 1)

	std::string file{"dataset_2.json"};
	//std::string file{};
	if(argc > 1) file = argv[1];

	//std::istringstream cmd("");
	std::string cmd;
	while (true) {
		allocation_optimize_NS::Problem *p{nullptr};

		if(!file.empty()) p = new allocation_optimize_NS::Problem(file);
		else if(argc <= 1) p = new allocation_optimize_NS::Problem();
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
			//std::pair<int, int> ans = allocation_optimize_NS::greedy(p);
			allocation_optimize_NS::Problem::Solution* ans = allocation_optimize_NS::greedy(p);
			timer.timerEnd(true);

			std::cout << "greedy optimization is end." << std::endl;

			delete ans;
			delete p;

			continue;

			//std::cin >> cmd;
			//return uint32_t(0);
		}
		else if(cmd == "b") {
			allocation_optimize_NS::Problem::Timer timer;
			std::cout << "Using brute force method:" << std::endl;
			std::list<allocation_optimize_NS::Problem::Solution> solutions = allocation_optimize_NS::brute_force(p);
			timer.timerEnd(true);
			for(auto s : solutions)  std::cout << s.getResultCost() << std::endl;
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
