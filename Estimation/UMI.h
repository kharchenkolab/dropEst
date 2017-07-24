#pragma once

#include <cstddef>
#include <vector>
#include <Tools/GtfRecord.h>

namespace Estimation
{
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
		private:
			char _mark;
		public:
			static const std::string DEFAULT_CODE;

			Mark(MarkType type = MarkType::NONE);

			void add(const Mark &mark);
			void add(MarkType type);
			void add(Tools::GtfRecord::RecordType type);
			bool check(MarkType type) const;
			bool match(const std::vector<Mark>) const;
			bool operator==(const MarkType &other) const;
			bool operator==(const Mark &other) const;

			static Mark get_by_code(char code);
			static std::vector<Mark> get_by_code(const std::string &code);
		};

	public:
		size_t read_count;
		Mark mark;

	public:
		UMI(size_t read_count = 0);
		void merge(const UMI& umi);
	};
}