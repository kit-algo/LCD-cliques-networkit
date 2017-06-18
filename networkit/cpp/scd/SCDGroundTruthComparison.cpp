#include "SCDGroundTruthComparison.h"
#include "../auxiliary/SignalHandling.h"
#include <string>

namespace NetworKit {
	SCDGroundTruthComparison::SCDGroundTruthComparison(const Graph &G, const Cover &groundTruth, const std::map<node, std::set<node>> &found, bool ignoreSeeds) : G(G), groundTruth(groundTruth), found(found), ignoreSeeds(ignoreSeeds) { };

	void SCDGroundTruthComparison::run() {
		Aux::SignalHandler handler;

		hasRun = false;

		jaccardScores.clear();
		f1Scores.clear();
		precisionScores.clear();
		recallScores.clear();

		std::vector<std::vector<node> > groundTruthSets (groundTruth.upperBound());
		std::vector<count> foundSizes(G.upperNodeIdBound(), 0), truthSizes(groundTruth.upperBound(), 0);

		G.forNodes([&](node u) {
			for (index s : groundTruth[u]) {
				groundTruthSets[s].emplace_back(u);
				++truthSizes[s];
			}
		});

		for (const auto &it : found) {
			for (node u : it.second) {
				if (G.hasNode(u)) {
					++foundSizes[it.first];
				}
			}
		}

		handler.assureRunning();

		std::map<index, count> overlap;

		for (const auto& found_it : found) {
			handler.assureRunning();

			node seed = found_it.first;

			if (!G.hasNode(seed)) {
				throw std::runtime_error("Error, the graph does not contain the seed node " + std::to_string(seed));
			}

			overlap.clear();

			const auto& foundNodes = found_it.second;

			auto allowedSubsets = groundTruth.subsetsOf(seed);

			for (node u : foundNodes) {
				if (G.hasNode(u)) {
					for (index s : groundTruth[u]) {
						if (ignoreSeeds || allowedSubsets.count(s) > 0) {
							++overlap[s];
						}
					}

				}
			}

			handler.assureRunning();

			double bestJaccard = 0;
			double bestF1 = 0;
			double bestPrecision = 0;
			double bestRecall = 0;

			for (auto o : overlap) {
				double currentj = o.second * 1.0 / (foundSizes[seed] + truthSizes[o.first] - overlap[o.first]);

				double recall = o.second * 1.0 / truthSizes[o.first];
				double precision = o.second * 1.0 / foundSizes[seed];

				double currentf1 = 0;
				if (precision > 0 && recall > 0) {
					currentf1 = 2 * (precision * recall) / (precision + recall);
				}

				if (currentj > bestJaccard) {
					bestJaccard = currentj;
				}

				if (currentf1 > bestF1) {
					bestF1 = currentf1;
				}

				if (recall > bestRecall) {
					bestRecall = recall;
				}

				if (precision > bestPrecision) {
					bestPrecision = precision;
				}

			}

			jaccardScores[seed] = bestJaccard;
			f1Scores[seed] = bestF1;
			precisionScores[seed] = bestPrecision;
			recallScores[seed] = bestRecall;
		}

		handler.assureRunning();

		auto avgScore = [](const std::map<index, double>& scores) -> double {
			double scoreSum = 0;
			for (auto a : scores) {
				scoreSum += a.second;
			}

			return scoreSum / scores.size();
		};


		averageJaccard = avgScore(jaccardScores);
		averageF1 = avgScore(f1Scores);
		averageRecall = avgScore(recallScores);
		averagePrecision = avgScore(precisionScores);

		handler.assureRunning();
		hasRun = true;
	};

	std::map< NetworKit::index, double > NetworKit::SCDGroundTruthComparison::getIndividualJaccard() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return jaccardScores;
	}

	std::map< NetworKit::index, double > NetworKit::SCDGroundTruthComparison::getIndividualF1() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return f1Scores;
	}

	std::map< NetworKit::index, double > NetworKit::SCDGroundTruthComparison::getIndividualRecall() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return recallScores;
	}

	std::map< NetworKit::index, double > NetworKit::SCDGroundTruthComparison::getIndividualPrecision() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return precisionScores;
	}

	double SCDGroundTruthComparison::getAverageJaccard() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return averageJaccard;
	}

	double SCDGroundTruthComparison::getAverageF1() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return averageF1;
	}

	double SCDGroundTruthComparison::getAverageRecall() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return averageRecall;
	}

	double SCDGroundTruthComparison::getAveragePrecision() {
		if (!hasRun) throw std::runtime_error("Error, run must be called first");
		return averagePrecision;
	}
};
