#ifndef MULTISORT_H
#define MULTISORT_H

#include <vector>
#include <algorithm>

namespace multiSort {
	template <typename T, typename Compare>
		std::vector<int> sort_permutation( std::vector<T> const& vec, Compare compare ) {
			std::vector<int> p(vec.size());
			std::iota(p.begin(), p.end(), 0);
			std::sort(p.begin(), p.end(),
					[&](int i, int j){ return compare(vec[i], vec[j]); });
			return p;
		}
	template <typename T>
		std::vector<int> sort_permutation( std::vector<T> const& vec ) {
			return sort_permutation( vec, std::less<T>() );
		}

	template <typename T>
		std::vector<T> apply_permutation( std::vector<T> const& vec, std::vector<int> const& p) {
			std::vector<T> sorted_vec(p.size());
			std::transform(p.begin(), p.end(), sorted_vec.begin(),
					[&](int i){ return vec[i]; });
			return sorted_vec;
		}
	template <typename T>
		void apply_permutation( std::vector<T>& vec, std::vector<int> const& p) {
			std::vector<T> sorted_vec(p.size());
			std::transform(p.begin(), p.end(), sorted_vec.begin(),
					[&](int i){ return vec[i]; });
			vec = std::move( sorted_vec );
		}

	template <typename T>
		void sort2( std::vector<T>& v1, std::vector<T>& v2 ) {
			std::vector<int> perm = sort_permutation( v1 );
			apply_permutation( v1, perm );
			apply_permutation( v2, perm );
		}
};

#endif
