// networkit-format
#include <unordered_map>

#include <networkit/auxiliary/UniformRandomSelector.hpp>
#include <networkit/clique/MaximalCliques.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/scd/CliqueDetect.hpp>

namespace NetworKit {

CliqueDetect::CliqueDetect(const Graph &G) : SelectiveCommunityDetector(G) {
    if (G.numberOfSelfLoops() > 0)
        throw std::runtime_error("CliqueDetect works only with simple graphs.");
}

std::set<node> CliqueDetect::expandOneCommunity(node s) {
    // 1. get the maximum clique in neighbors(s) + s
    std::set<node> result;
    result.insert(s);

    std::vector<node> sn(G->neighborRange(s).begin(), G->neighborRange(s).end());
    std::vector<edgeweight> neighborWeights;
    if (G->isWeighted()) {
        neighborWeights.reserve(sn.size());

        for (auto it : G->weightNeighborRange(s)) {
            neighborWeights.push_back(it.second);
        }
    }

    if (!sn.empty()) {
        for (node u : getMaximumWeightClique(sn, neighborWeights)) {
            result.insert(sn[u]);
        }
    }

    return result;
}

std::set<node> CliqueDetect::expandOneCommunity(const std::set<node> &seeds) {
    std::set<node> result(seeds);

    if (!seeds.empty()) {
        // Find neighbors of the seeds that are neighbors of all seed nodes
        // First, candidates are neighbors of the first seed node
        std::unordered_map<node, std::pair<count, edgeweight>> sn;
        for (auto neighborIt : G->weightNeighborRange(*seeds.begin())) {
            node v = neighborIt.first;
            edgeweight weight = neighborIt.second;
            if (seeds.find(v) == seeds.end()) {
                sn.insert({v, {1u, weight}});
            }
        }

        // Count for each neighbor to how many other seed nodes it is adjacent
        for (auto it = ++seeds.begin(); it != seeds.end(); ++it) {
            for (auto neighborIt : G->weightNeighborRange(*it)) {
                node v = neighborIt.first;
                edgeweight weight = neighborIt.second;
                auto it = sn.find(v);
                if (it != sn.end()) {
                    ++it->second.first;
                    it->second.second += weight;
                }
            }
        }

        // Take those neighbors that are adjacent to all seed nodes
        std::vector<node> neighbors;
        std::vector<edgeweight> neighborWeights;
        neighbors.reserve(sn.size());
        if (G->isWeighted()) {
            neighborWeights.reserve(sn.size());
        }

        for (auto it : sn) {
            assert(it.second.first <= seeds.size() && "Graph has probably multi-edges!");
            if (it.second.first == seeds.size()) {
                neighbors.push_back(it.first);
                if (G->isWeighted()) {
                    neighborWeights.push_back(it.second.second);
                }
            }
        }

        if (!neighbors.empty()) {
            const std::vector<node> clique = getMaximumWeightClique(neighbors, neighborWeights);

            for (node u : clique) {
                result.insert(neighbors[u]);
            }
        }
    }

    return result;
}

std::vector<node>
CliqueDetect::getMaximumWeightClique(const std::vector<node> &nodes,
                                     const std::vector<edgeweight> &seedToNodeWeight) const {

    Graph S = createSubgraphFromNodes(nodes);
    std::vector<node> maxClique;

    Aux::UniformRandomSelector selector;
    if (!G->isWeighted()) {
        // Select a maximum clique uniformly at random
        MaximalCliques mc(S, [&](const std::vector<node> &clique) {
            if (clique.size() < maxClique.size())
                return;
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
        // Select a maximum weight clique uniformly at random
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
                for (auto neighborIt : S.weightNeighborRange(u)) {
                    cliqueWeight += inClique[neighborIt.first] * neighborIt.second;
                }

                // Add the weight from u to the seed node(s)
                cliqueWeight += seedToNodeWeight[u];

                // Reset inClique[u] to avoid counting edges twice
                inClique[u] = false;
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
        });

        mc.run();
    }

    return maxClique;
}

Graph CliqueDetect::createSubgraphFromNodes(const std::vector<node> &nodes) const {
    GraphBuilder gbuilder(nodes.size(), G->isWeighted());
    std::unordered_map<node, node> reverseMapping;

    node i = 0;
    for (node u : nodes) {
        reverseMapping[u] = i;
        i += 1;
    }
    for (node u : nodes) {
        node lu = reverseMapping[u];
        for (auto neighborIt : G->weightNeighborRange(u)) {
            auto lv = reverseMapping.find(neighborIt.first);
            if (lv != reverseMapping.end())
                gbuilder.addHalfEdge(lu, lv->second, neighborIt.second);
        }
    }

    return gbuilder.toGraph(false);
}

} /* namespace NetworKit */
