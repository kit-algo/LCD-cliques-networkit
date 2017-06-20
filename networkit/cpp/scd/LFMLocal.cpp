#include "LFMLocal.h"

#include <cmath>
#include <unordered_map>
#include "../structures/LocalCommunity.h"
#include "../auxiliary/UniformRandomSelector.h"

namespace NetworKit {

LFMLocal::LFMLocal(const Graph& G, const double alpha): SelectiveCommunityDetector(G), alpha(alpha) {
}


std::set<node> LFMLocal::expandOneCommunity(const std::set<node>& seeds) {
	LocalCommunity<true, false, true> community(G);
	using shell_info_t = LocalCommunity<true, false, true>::ShellInfo;
	using community_info_t = LocalCommunity<true, false, true>::CommunityInfo;

	for (node u : seeds) {
		community.addNode(u);
	}

	if (community.internalEdgeWeight() + community.cut() == 0) return community.toSet();

	auto quality = [this](double internalEdgeWeight, double cut) {
		return (2 * internalEdgeWeight) / std::pow(static_cast<double>(2 * internalEdgeWeight + cut), alpha);
	};

	node uMax;
	double qMax, currentQ = quality(community.internalEdgeWeight(), community.cut());
	do {
		assert(std::abs(currentQ - quality(community.internalEdgeWeight(), community.cut())) < 0.00001);
		uMax = none;
		qMax = 0.0;

		// find node with max score

		Aux::UniformRandomSelector selector;
		community.forShellNodes([&](node u, const shell_info_t& u_info) {
			double uQ = quality(community.internalEdgeWeight() + u_info.intDeg.get(),
					community.cut() - u_info.intDeg.get() + u_info.extDeg.get()) - currentQ;

			assert(!std::isnan(uQ) && !std::isinf(uQ));
			TRACE(u, " diff: ", uQ);
			if (uQ > qMax) {
				qMax = uQ;
				uMax = u;
				selector.reset();
			} else if (uQ == qMax && selector.add_element()) {
				uMax = u;
			}
		});

		if (uMax != none) {
			community.addNode(uMax);
			currentQ += qMax;

			// delete nodes from result with neg. score (except node s)
			Aux::UniformRandomSelector removalSelector;
			node uMin;
			do {
				uMin = none;
				double qMin = 0;

				assert(std::abs(currentQ - quality(community.internalEdgeWeight(), community.cut())) < 0.00001);

				community.forCommunityNodes([&](node u, const community_info_t& u_info) {
					double uQ = currentQ - quality(community.internalEdgeWeight() - u_info.intDeg.get(), community.cut() + u_info.intDeg.get() - u_info.extDeg.get());

					TRACE(u, " diff removal: ", uQ);

					if (uQ < qMin && seeds.find(u) == seeds.end()) {
						uMin = u;
						qMin = uQ;
						removalSelector.reset();
					} else if (uQ == qMin && seeds.find(u) == seeds.end() && removalSelector.add_element()) {
						uMin = u;
					}
				});

				if (uMin != none) {
					currentQ -= qMin;
					community.removeNode(uMin);
				}
			} while (uMin != none);
		}

	} while (uMax != none);

	return community.toSet();
}

} // namespace NetworKit
