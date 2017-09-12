#pragma once

#include <cstddef>
#include <queue>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_lock_t int
	#define omp_set_lock(x)
	#define omp_unset_lock(x)
	#define omp_init_lock(x)
	#define omp_destroy_lock(x)
#endif

#include "Logs.h"

namespace Tools
{
	template <typename T>
	class OmpConcurrentQueue
	{
	private:
		class RaiiOmpLock
		{
		private:
			omp_lock_t *_lock;

		public:
			RaiiOmpLock(omp_lock_t *lock)
			{
				omp_set_lock(lock);
				this->_lock = lock;
			}

			~RaiiOmpLock()
			{
				omp_unset_lock(this->_lock);
			}
		};

	private:
		std::queue<T> _queue;
		omp_lock_t _lock;

		size_t _size;
		const size_t _max_size;

	public:
		OmpConcurrentQueue(size_t max_size)
			: _size(0)
			, _max_size(max_size)
		{
			omp_init_lock(&this->_lock);
		}

		~OmpConcurrentQueue()
		{
			omp_destroy_lock(&this->_lock);
		}

		bool pop(T& item)
		{
			RaiiOmpLock l(&this->_lock);

			if (this->_size == 0)
				return false;

			item = this->_queue.front();
			this->_queue.pop();
			this->_size--;

			return true;
		}

		void push(const T& item)
		{
			RaiiOmpLock l(&this->_lock);

			this->_queue.push(item);
			this->_size++;
		}

		void push(T&& item)
		{
			RaiiOmpLock l(&this->_lock);

			this->_queue.push(std::move(item));
			this->_size++;
		}

		bool full() const
		{
			size_t size;
			#pragma omp atomic read
			size = this->_size;
			return size >= this->_max_size;
		}

		bool empty() const
		{
			size_t size;
			#pragma omp atomic read
			size = this->_size;
			return size == 0;
		}
	};
}

