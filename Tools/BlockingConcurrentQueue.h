#pragma once

#include <cstddef>
#include <queue>
#include <mutex>
#include <atomic>

#include "Logs.h"

namespace Tools
{
	template <typename T>
	class BlockingConcurrentQueue
	{
	private:
		typedef std::mutex mutex_t;
		typedef std::lock_guard<mutex_t> lock_t;

	private:
		std::queue<T> _queue;
		std::mutex _lock;

		std::atomic<size_t> _size;
		const size_t _max_size;

	public:
		BlockingConcurrentQueue(size_t max_size)
			: _queue()
			, _size(0)
			, _max_size(max_size)
		{}

		bool pop(T& item)
		{
			lock_t l(this->_lock);

			if (this->_size == 0)
				return false;

			item = this->_queue.front();
			this->_queue.pop();
			this->_size--;

			return true;
		}

		void push(const T& item)
		{
			lock_t l(this->_lock);

			this->_queue.push(item);
			this->_size++;
		}

		void push(T&& item)
		{
			lock_t l(this->_lock);

			this->_queue.push(std::move(item));
			this->_size++;
		}

		bool full() const // shouldn't block the queue. It's just a suggested size.
		{
			return this->_size >= this->_max_size;
		}

		bool empty() const
		{
			return this->_size == 0;
		}
	};
}

