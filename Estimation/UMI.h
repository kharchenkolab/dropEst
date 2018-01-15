#pragma once

#include <cstddef>
#include <vector>
#include <Tools/GeneAnnotation/GtfRecord.h>

namespace Estimation
{
	class ReadInfo;
	class UMI
	{
	public:
		class Mark
		{
		public:
			enum MarkType
			{
				NONE = 0,
				HAS_NOT_ANNOTATED = 1,
				HAS_EXONS = 2,
				HAS_INTRONS = 4
			};

			using query_t = std::vector<Mark>;

		private:
			char _mark;

		public:
			static const std::string DEFAULT_CODE;

			explicit Mark(MarkType type = MarkType::NONE);

			void add(const Mark &mark);
			void add(MarkType type);
			void add(Tools::GeneAnnotation::GtfRecord::RecordType type);
			bool check(MarkType type) const;
			bool match(const std::vector<Mark> &match_levels) const;
			bool operator==(const MarkType &other) const;
			bool operator==(const Mark &other) const;

			static Mark get_by_code(char code);
			static std::vector<Mark> get_by_code(const std::string &code);
		};

		using quality_t = std::vector<unsigned>;

	private:
		size_t _read_count;
		Mark _mark;
		quality_t _sum_quality;

	public:
		explicit UMI(size_t quality_length = 0, size_t read_count = 0);

		size_t read_count() const;
		const Mark& mark() const;
		std::vector<double> mean_quality() const;

		void merge(const UMI& umi);
		void add_read(const ReadInfo &read_info);
	};
}
