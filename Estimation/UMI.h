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

			typedef std::vector<Mark> query_t;

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

	private:
		size_t _read_count;
		Mark _mark;

	public:
		UMI(size_t read_count = 0);

		size_t read_count() const;
		const Mark& mark() const;

		void merge(const UMI& umi);
		void add_read(Mark mark);
	};
}
