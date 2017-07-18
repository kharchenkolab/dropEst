#pragma once

#include <string>
#include <vector>

namespace Tools
{
	class IndexedValue;
}

namespace TestEstimator
{
	struct testFillDistances;
	struct testRealNeighboursCbs;
	struct testRealNeighbours;
	struct testBarcodesFile;
	struct testConstLengthBarcodesFile;
}

namespace Estimation
{
namespace Merge
{
namespace BarcodesParsing
{
	class BarcodesParser
	{
		friend struct TestEstimator::testFillDistances;
		friend struct TestEstimator::testRealNeighboursCbs;
		friend struct TestEstimator::testRealNeighbours;
		friend struct TestEstimator::testBarcodesFile;
		friend struct TestEstimator::testConstLengthBarcodesFile;

	public:
		struct BarcodesDistance
		{
			std::vector<size_t> barcode_part_inds;
			unsigned edit_distance;

			BarcodesDistance(const std::vector<size_t> &barcodes_inds, unsigned edit_distance);
		};

	protected:
		typedef std::vector<std::string> barcodes_list_t;
		typedef std::vector<barcodes_list_t> barcode_parts_list_t;
		typedef std::vector<std::vector<Tools::IndexedValue>> edit_distance_parts_list_t;
		typedef std::vector<BarcodesDistance> barcodes_distance_list_t;

	private:
		const std::string _barcodes_filename;
		barcode_parts_list_t _barcodes;

	protected:
		static const int MAX_REAL_MERGE_EDIT_DISTANCE=5;
		const std::string& barcode(size_t part_ind, size_t barcode_ind) const;
		const size_t barcode_parts_num() const;

	private:
		edit_distance_parts_list_t get_distances_to_barcode(const std::string &barcode) const;
		void push_remaining_dists(edit_distance_parts_list_t::const_iterator begin,
								  edit_distance_parts_list_t::const_iterator end, unsigned edit_distance,
								  const std::vector<size_t> &barcodes_inds, barcodes_distance_list_t &res) const;

	protected:
		virtual barcodes_list_t split_barcode(const std::string &barcode) const = 0;
		virtual barcode_parts_list_t get_barcodes_list(const std::string &barcodes_filename) const = 0;

	public:
		BarcodesParser(const std::string &barcodes_filename);
		virtual void init();
		virtual void release();

		virtual std::string get_barcode(const std::vector<size_t> &barcode_part_inds) const;
		virtual barcodes_distance_list_t get_real_neighbour_cbs(const std::string &barcode) const;
	};
}
}
}