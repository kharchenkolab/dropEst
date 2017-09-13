#pragma once

#include <cstddef>
#include <queue>
#include <mutex>
#include <atomic>

#include "Logs.h"
#include "Tools/readerwriterqueue.h"

namespace Tools
{
	template <typename T>
	class ScSpConcurrentQueue
	{
	private:
		moodycamel::ReaderWriterQueue<T> _queue;

		std::atomic<size_t> _size;
		const size_t _max_size;

	public:
		ScSpConcurrentQueue(size_t max_size)
			: _queue()
			, _size(0)
			, _max_size(max_size)
		{}

		bool pop(T& item)
		{
			if (!this->_queue.try_dequeue(item))
				return false;

			this->_size--;

			return true;
		}

		void push(const T& item)
		{
			this->_queue.enqueue(item);
			this->_size++;
		}

		void push(T&& item)
		{
			this->_queue.enqueue(std::move(item));
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

