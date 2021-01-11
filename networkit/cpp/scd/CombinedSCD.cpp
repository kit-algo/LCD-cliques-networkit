// networkit-format
#include <networkit/scd/CombinedSCD.hpp>

namespace NetworKit {

CombinedSCD::CombinedSCD(const Graph &G, SelectiveCommunityDetector &first,
                         SelectiveCommunityDetector &second)
    : SelectiveCommunityDetector(G), first(first), second(second) {}

std::set<node> CombinedSCD::expandOneCommunity(node s) {
    return second.expandOneCommunity(first.expandOneCommunity(s));
}

std::set<node> CombinedSCD::expandOneCommunity(const std::set<node> &s) {
    return second.expandOneCommunity(first.expandOneCommunity(s));
}

} // namespace NetworKit
