#ifndef MINMAX_HEAP2_H_
#define MINMAX_HEAP2_H_

#ifndef NDEBUG
#define INTERNAL_HEAP_CHECKS
#endif

#include <utility>
#include <vector>
#include <cassert>
#include <functional>
#include <algorithm>
#include <unordered_map>

/**
 * Minimum/maximum heap with dynamic expansion of the size of the heap.
 *
 * @author Ben Strasser <strasser@kit.edu>
 * @author Michael Hamann <michael.hamann@kit.edu>
 */
namespace Aux {

template<class keyT, unsigned int heapArity = 4, class keyOrderT = std::less<keyT> >
class MinHeap{
public:
	typedef keyT key_type;
	typedef keyOrderT key_order_type;

	explicit MinHeap(unsigned int capacity, key_order_type order = key_order_type()):
		heap_end(0), heap(), id_pos(), order(std::move(order)), contained_flags() {
			
		heap.reserve(capacity);

		check_id_invariants();
		check_order_invariants();
	}

	explicit MinHeap(key_order_type order = key_order_type()):
		heap_end(0), heap(), id_pos(), order(std::move(order)), contained_flags() {
		
		check_id_invariants();
		check_order_invariants();
	}

	//! Takes a range of std::pair<int, key_type> elements (or any other struct with first and second) 
	//! and builds a heap containing these elements.
	template<class Range>
	void fill(const Range& range) {
		clear();

		for (auto r : range) {
			if (r >= id_pos.size()) {
				id_pos.resize(r + 1);
				contained_flags.resize(r + 1, false);
			}
			assert(0 <= r.first && r.first < (int)id_pos.size() && "range must contain valid id");
			heap.resize(heap_end + 1);

			heap[heap_end].id = r.first;
			heap[heap_end].key = r.second;

			id_pos[r.first] = heap_end;
			contained_flags[r.first] = true;
			heap_end += 1;
		}

		rebuildHeap();

		check_id_invariants();
		check_order_invariants();
	}


	void clear() {
		contained_flags.clear();
		id_pos.clear();
		heap.clear();
		heap_end = 0;

		check_id_invariants();
		check_order_invariants();
	}

	void reorder(key_order_type new_order) {
		check_id_invariants();
		check_order_invariants();

		order = std::move(new_order);
		rebuildHeap();

		check_id_invariants();
		check_order_invariants();
	}

	void reset(key_order_type new_order) {
		clear();
		order = std::move(new_order);
	}

	bool empty() const {
		assert(heap.size() == heap_end);

		check_id_invariants();
		check_order_invariants();

		return heap_end == 0;
	}

	size_t size() const {
		assert(heap.size() == heap_end);

		return heap_end;
	}

	bool contains(unsigned int id) const {
		if (id >= contained_flags.size()) return false;

		check_id_invariants();
		check_order_invariants();
		
		return contained_flags[id];
	}

	const key_type& getKey(unsigned int id) const {
		assert(contains(id) && "id is contained");

		check_id_invariants();
		check_order_invariants();

		return heap[id_pos.at(id)].key;
	}

	void push(unsigned int id, key_type key) {
		if (id >= id_pos.size()) {
			id_pos.resize(id + 1);
			contained_flags.resize(id + 1, false);
		}
		assert(0 <= id && id < id_pos.size() && "id is in range");
		assert(!contains(id) && "can not push an already existing id");
		assert(heap.size() == heap_end);

		contained_flags[id] = true;

		size_t new_pos = heap_end;
		heap_end += 1;
		heap.emplace_back(IdKeyPair {id, std::move(key)});
		assert(heap.size() == heap_end);

		id_pos[id] = new_pos;
		moveUp(new_pos);

		check_id_invariants();
		check_order_invariants();
	}

	bool pushOrDecreaseKey(unsigned int id, key_type key) {
		check_id_invariants();
		check_order_invariants();
		
		if (!contains(id)) {
			push(id, key);
			return true;
		} else if (order(key, heap[id_pos[id]].key)) {
				heap[id_pos[id]].key = std::move(key);
				moveUp(id_pos[id]);

				check_id_invariants();
				check_order_invariants();
		
				return true;
		}
		return false;
	}

	bool pushOrIncreaseKey(unsigned int id, key_type key) {
		check_id_invariants();
		check_order_invariants();
		
		if (!contains(id)) {
			push(id, key);
			return true;
		} else if (order(heap[id_pos[id]].key, key)) {
				heap[id_pos[id]].key = std::move(key);
				moveDown(id_pos[id]);
				
				check_id_invariants();
				check_order_invariants();
		
				return true;
		}
		return false;
	}

	bool pushOrSetKey(unsigned int id, key_type key) {
		check_id_invariants();
		check_order_invariants();
		
		if (!contains(id)) {
			push(id, key);
			return true;
		} else if (order(heap[id_pos[id]].key, key)) {
			heap[id_pos[id]].key = std::move(key);
			moveDown(id_pos[id]);
				
			check_id_invariants();
			check_order_invariants();
		
			return true;
		} else if (order(key, heap[id_pos[id]].key)) {
			heap[id_pos[id]].key = std::move(key);
			moveUp(id_pos[id]);

			check_id_invariants();
			check_order_invariants();
		
			return true;
		}
		return false;
	}

	key_type peekMinKey() const {
		assert(!empty() && "heap is not empty");
		
		check_id_invariants();
		check_order_invariants();
		
		return heap[0].key;
	}

	unsigned int peekMinId() const {
		assert(!empty() && "heap is not empty");

		check_id_invariants();
		check_order_invariants();
		
		return heap[0].id;
	}

	unsigned int pop() {
		assert(!empty() && "heap is not empty");
		assert(heap.size() == heap_end);

		check_id_invariants();
		check_order_invariants();
		
		if (heap_end == 1) {
			heap_end = 0;
			contained_flags[heap[0].id] = false;
			heap.resize(heap_end);

			check_id_invariants();
			check_order_invariants();
			
			return heap[0].id;
		} else {	
			unsigned int ret = heap[0].id;
			heap_end -= 1;

			#ifdef INTERNAL_HEAP_CHECKS
			assert(ret != heap[heap_end].id);
			#endif
			
			heap[0].id = heap[heap_end].id;
			heap[0].key = std::move(heap[heap_end].key);
			id_pos[heap[0].id] = 0;
			moveDown(0);
			contained_flags[ret] = false;
			heap.resize(heap_end);
			
			check_id_invariants();
			check_order_invariants();

			return ret;
		}
	}

private:

	static unsigned int parent(unsigned int pos) {
		assert(pos != 0);
		return (pos - 1) / heapArity;
	}

	static unsigned int childrenBegin(unsigned int pos) {
		return heapArity * pos + 1;
	}

	static unsigned int childrenEnd(unsigned int pos) {
		return heapArity * (pos + 1) + 1;
	}

	void rebuildHeap() {
		for (int i = heap_end - 1; i >= 0; --i)
			moveDown(i);
	}

	void moveUp(unsigned int pos) {
		if (pos != 0) {
			key_type key = std::move(heap[pos].key);
			unsigned int id = heap[pos].id;

			unsigned int parent_pos = parent(pos);
			while(order(key, heap[parent_pos].key)) {
				heap[pos].id = heap[parent_pos].id;
				heap[pos].key = std::move(heap[parent_pos].key);
				id_pos[heap[parent_pos].id] = pos;

				pos = parent_pos;
				if (pos == 0)
					break;
				parent_pos = parent(pos);
			}

			heap[pos].id = id;
			heap[pos].key = std::move(key);
			id_pos[id] = pos;
		}
	}

	void moveDown(unsigned int pos) {
		key_type key = std::move(heap[pos].key);
		unsigned int id = heap[pos].id;

		for (;;) {
			unsigned int begin = std::min(static_cast<unsigned int>(heap_end), childrenBegin(pos));
			unsigned int end = std::min(static_cast<unsigned int>(heap_end), childrenEnd(pos));

			if (begin == end)
				break;

			unsigned int min_child_pos = begin;
			for (unsigned int i = begin + 1; i < end; ++i) {
				if (order(heap[i].key, heap[min_child_pos].key))
					min_child_pos = i;
			}

			if (!order(heap[min_child_pos].key, key))
				break;

			heap[pos].id = heap[min_child_pos].id;
			heap[pos].key = std::move(heap[min_child_pos].key);
			id_pos[heap[min_child_pos].id] = pos;

			pos = min_child_pos;
		}
		heap[pos].id = id;
		heap[pos].key = std::move(key);
		id_pos[id] = pos;
	}

	struct IdKeyPair {
		unsigned int id;
		key_type key;
	};

	size_t heap_end;
	std::vector<IdKeyPair> heap;
	std::vector<size_t> id_pos;
	key_order_type order;
	std::vector<bool> contained_flags;

	void check_id_invariants() const {
		#ifdef INTERNAL_HEAP_CHECKS
		for (unsigned int i = 0; i < heap_end; ++i) {
			//assert(heap[i].id != -1);
			assert(0 <= heap[i].id);
			//assert(heap[i].id < id_pos.size());
			assert(contained_flags.at(heap[i].id));
			assert(id_pos.at(heap[i].id) == i);
		}

		for (unsigned int i = 0; i < id_pos.size(); ++i) {
			if (contained_flags[i]) {
				assert(0 <= id_pos[i]);
				assert(id_pos[i] < heap_end);
				assert(heap[id_pos[i]].id == i);
			}
		}

		//assert(heap.size() == id_pos.size());
		//assert(heap.size() == contained_flags.size());
		#endif
	}

	void check_order_invariants() const {
		#ifdef INTERNAL_HEAP_CHECKS
		// TODO fix
		//for (unsigned int i = 1; i < heap_end; ++i)
		//	assert(!order(heap[i].key, heap[parent(i)].key));
		#endif
	}
};

/**
 * MaxHeap is a 'heapArity'-ary Heap with items (int, 'keyT'), which is ordered
 * by 'keyT' with an order of type 'keyOrderT'.
 *
 * This class is based on MinHeap with inverse order. For further comments
 * @see MinHeap.
 */
template<class keyT, unsigned int heapArity = 4, class keyOrderT = std::less<keyT>>
class MaxHeap{
public:
	typedef keyT key_type;
	typedef keyOrderT key_order_type;

	explicit MaxHeap(int id_count, key_order_type order = key_order_type())
		:heap(id_count, inverseOrder(std::move(order))) {}

	explicit MaxHeap(key_order_type order = key_order_type())
		:heap(inverseOrder(std::move(order))) {}

	void clear() {
		heap.clear();
	}

	template<class Range>
	void fill(const Range& range) {
		heap.fill(range);
	}
	
	void reorder(key_order_type new_order) {
		heap.reorder(inverseOrder(std::move(new_order)));
	}

/*	void reset(int new_id_count = 0, key_order_type new_order = key_order_type()) {
		heap.reset(new_id_count, inverseOrder(std::move(new_order)));
	} */

	void reset(key_order_type new_order) {
		heap.reset(std::move(new_order));
	}

	bool empty() const {
		return heap.empty();
	}

	size_t size() const {
		return heap.size();
	}

	bool contains(unsigned int id) const {
		return heap.contains(id);
	}

	const key_type& getKey(unsigned int id) const {
		return heap.getKey(id);
	}

	void push(unsigned int id, key_type key) {
		heap.push(id, std::move(key));
	}

	bool pushOrDecreaseKey(unsigned int id, key_type key) {
		return heap.pushOrIncreaseKey(id, std::move(key));
	}

	bool pushOrIncreaseKey(unsigned int id, key_type key) {
		return heap.pushOrDecreaseKey(id, std::move(key));
	}

	bool pushOrSetKey(unsigned int id, key_type key) {
		return heap.pushOrSetKey(id, std::move(key));
	}

	key_type peekMaxKey() const {
		return heap.peekMinKey();
	}

	unsigned int peekMaxId() const {
		return heap.peekMinId();
	}

	unsigned int pop() {
		return heap.pop();
	}

private:
	struct inverseOrder{
		inverseOrder() {}

		inverseOrder(key_order_type order):order(std::move(order)) {}
		
		bool operator()(const key_type&l, const key_type&r) {
			return order(r, l);
		}
		
		key_order_type order;
	};
	
	MinHeap<key_type, heapArity, inverseOrder> heap;
};

} /* namespace Aux */

#endif /* MINMAX_HEAP2_H_ */
