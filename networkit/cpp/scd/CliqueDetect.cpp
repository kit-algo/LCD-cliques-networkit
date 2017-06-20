#include "CliqueDetect.h"
#include "../clique/MaximalCliques.h"
#include "../graph/GraphBuilder.h"
#include "../auxiliary/UniformRandomSelector.h"
#include <unordered_map>

namespace NetworKit {

CliqueDetect::CliqueDetect(const Graph& G) : SelectiveCommunityDetector(G) {
	if (G.numberOfSelfLoops() > 0) throw std::runtime_error("CliqueDetect works only with simple graphs.");
}

std::set<node> CliqueDetect::expandOneCommunity(node s) {
	// 1. get the maximum clique in neighbors(s) + s
	std::set<node> result;
	result.insert(s);

	const std::vector<node> sn = G.neighbors(s);

	if (!sn.empty()) {
		for (node u : getMaximumWeightClique(sn, [s](node v) {return v == s;})) {
			result.insert(sn[u]);
		}
	}

	return result;
}

std::set<node> CliqueDetect::expandOneCommunity(const std::set<node>& seeds) {
	std::set<node> result(seeds);

	if (!seeds.empty()) {
		std::unordered_map<node, count> sn;
		G.forNeighborsOf(*seeds.begin(), [&](node v) {
			if (seeds.find(v) == seeds.end()) {
				sn.insert({v, 1u});
			}
		});

		for (auto it = ++seeds.begin(); it != seeds.end(); ++it) {
			G.forNeighborsOf(*it, [&](node v) {
				auto it = sn.find(v);
				if (it != sn.end()) {
					++it->second;
				}
			});
		}

		std::vector<node> neighbors;
		neighbors.reserve(sn.size());
		for (auto it : sn) {
			assert(it.second <= seeds.size() && "Graph has probably multi-edges!");
			if (it.second == seeds.size()) {
				neighbors.push_back(it.first);
			}
		}

		if (neighbors.size() > 0) {
			const std::vector<node> clique = getMaximumWeightClique(neighbors, [&](node v) { return seeds.find(v) != seeds.end(); });

			for (node u : clique) {
				result.insert(neighbors[u]);
			}
		}
	}

	return result;
}

template <typename F>
std::vector<node> CliqueDetect::getMaximumWeightClique(const std::vector<node>& nodes, F is_seed) const {
	Graph S = createSubgraphFromNodes(nodes);
	std::vector<node> maxClique;

	Aux::UniformRandomSelector selector;
	if (!G.isWeighted()) {
		MaximalCliques mc(S, [&](const std::vector<node>& clique) {
			if (clique.size() < maxClique.size()) return;
			if (clique.size() > maxClique.size()) {
				maxClique = clique;
				selector.reset();
			} else {
				if (selector.add_element()) {
					maxClique = clique;
				}
			}
		});
		mc.run();
	} else {
		double maxWeight = 0.0;
		assert(S.numberOfNodes() == S.upperNodeIdBound());
		std::vector<bool> inClique(S.upperNodeIdBound(), false);

		// clique with the greatest sum of its edgeweights
		MaximalCliques mc(S, [&](const std::vector<node> &clique) {
			double cliqueWeight = 0.0;

			for (node u : clique) {
				inClique[u] = true;
			}

			for (node u : clique) {
				S.forNeighborsOf(u, [&] (node, node v, edgeweight weight) {
					cliqueWeight += inClique[v] * weight/2;
					cliqueWeight += is_seed(v) * weight;
				});
			}

			if (cliqueWeight > maxWeight) {
				maxWeight = cliqueWeight;
				maxClique = clique;
				selector.reset();
			} else if (cliqueWeight == maxWeight) {
				if (selector.add_element()) {
					maxClique = clique;
				}
			}

			for (node u : clique) {
				inClique[u] = false;
			}
		});

		mc.run();
	}

	return maxClique;
}


Graph CliqueDetect::createSubgraphFromNodes(const std::vector<node>& nodes) const {
	GraphBuilder gbuilder(nodes.size(), true);
	std::unordered_map<node, node> reverseMapping;

	node i = 0;
	for (node u : nodes) {
		reverseMapping[u] = i;
		i += 1;
	}
	for (node u : nodes) {
		node lu = reverseMapping[u];
		G.forNeighborsOf(u, [&] (node, node v, edgeweight ew) {
			auto lv = reverseMapping.find(v);
			if (lv != reverseMapping.end())
				gbuilder.addHalfEdge(lu, lv->second, ew);
		});
	}

	return gbuilder.toGraph(false);
}

} /* namespace NetworKit */
