#include <iterator>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/scd/RandomBFS.hpp>

namespace NetworKit {

RandomBFS::RandomBFS(const Graph &G, const Cover &C)
    : SelectiveCommunityDetector(G), C(&C), subsetSizes(C.subsetSizeMap()) {}

std::set<node> RandomBFS::expandOneCommunity(const std::set<node> &s) {
    // allow as many nodes in the community as seeds are given
    count com_size = s.size();

    { // select a random community of s and get its size
        std::set<index> gs(C->subsetsOf(*s.begin()));

        for (auto it = ++s.begin(); it != s.end(); ++it) {
            const std::set<index> additional_coms(C->subsetsOf(*it));
            for (auto it = gs.begin(); it != gs.end();) {
                if (additional_coms.find(*it) == additional_coms.end()) {
                    it = gs.erase(it);
                } else {
                    ++it;
                }
            }
        }

        if (!gs.empty()) {
            index i = Aux::Random::index(gs.size());
            auto it = gs.begin();
            std::advance(it, i);
            com_size = subsetSizes[*it];
        }
    }

    std::set<node> result;

    std::vector<node> current_level, next_level;

    for (node u : s) {
        current_level.push_back(u);
    }

    while (result.size() < com_size && !current_level.empty()) {
        next_level.clear();

        if (current_level.size() < result.size() + com_size) {
            for (node u : current_level) {
                result.insert(u);
            }
        } else {
            std::shuffle(current_level.begin(), current_level.end(), Aux::Random::getURNG());

            for (auto it = current_level.begin(); result.size() < com_size; ++it) {
                assert(it != current_level.end());
                result.insert(*it);
            }

            break;
        }

        for (node u : current_level) {
            G->forNeighborsOf(u, [&](node v) {
                if (result.count(v) == 0) {
                    next_level.push_back(v);
                }
            });
        }

        std::sort(next_level.begin(), next_level.end());
        auto last = std::unique(next_level.begin(), next_level.end());
        next_level.erase(last, next_level.end());

        std::swap(current_level, next_level);
    }

    return result;
}

} // namespace NetworKit
