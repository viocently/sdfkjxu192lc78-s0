// Copyright(c) 2020 Robert Vaser
// Combination of ThreadPool implementation by progschj and
//   task stealing by Sean Parent

#ifndef THREAD_POOL_THREAD_POOL_HPP_
#define THREAD_POOL_THREAD_POOL_HPP_

#include <iostream>
#include <algorithm>
#include <atomic>
#include <functional>
#include <future>  // NOLINT
#include <memory>
#include <queue>
#include <thread>  // NOLINT
#include <unordered_map>
#include <utility>
#include <vector>

namespace thread_pool {

    class ThreadPool {
    public:
        explicit ThreadPool(
            std::size_t num_threads = std::thread::hardware_concurrency())
            : threads_(),
            thread_map_(),
            queues_(1UL > num_threads ? 1UL : num_threads),
            task_id_(0) {
            for (std::size_t i = 0; i != queues_.size(); ++i) {
                threads_.emplace_back([this, i]() -> void { Task(i); });
                thread_map_.emplace(threads_.back().get_id(), i);
            }
        }

        ThreadPool(const ThreadPool&) = delete;
        ThreadPool& operator=(const ThreadPool&) = delete;

        ThreadPool(ThreadPool&&) = delete;
        ThreadPool& operator=(ThreadPool&&) = delete;

        ~ThreadPool() {
            for (auto& it : queues_) {
                it.Done();
            }
            for (auto& it : threads_) {
                it.join();
            }
        }

        std::size_t num_threads() const {
            return threads_.size();
        }

        const std::unordered_map<std::thread::id, std::size_t>& thread_map() const {
            return thread_map_;
        }

        template<typename T, typename... Ts>
        auto Submit(T&& routine, Ts&&... params)
            -> std::future<typename std::result_of<T(Ts...)>::type> {
            auto task = std::make_shared<std::packaged_task<typename std::result_of<T(Ts...)>::type()>>(  // NOLINT
                std::bind(std::forward<T>(routine), std::forward<Ts>(params)...));
            auto task_result = task->get_future();
            auto task_wrapper = [task]() {
                (*task)();
            };

            auto task_id = task_id_++;
            bool is_submitted = false;
            for (std::size_t i = 0; i != queues_.size() * 42; ++i) {
                if (queues_[(task_id + i) % queues_.size()].TryPush(task_wrapper)) {
                    is_submitted = true;
                    break;
                }
            }
            if (!is_submitted) {
                queues_[task_id % queues_.size()].Push(task_wrapper);
            }

            return task_result;
        }

    private:
        void Task(std::size_t thread_id) {
            while (true) {
                std::function<void()> task;

                for (std::size_t i = 0; i != queues_.size(); ++i) {
                    if (queues_[(thread_id + i) % queues_.size()].TryPop(&task)) {
                        break;
                    }
                }
                if (!task && !queues_[thread_id].Pop(&task)) {
                    break;
                }

                task();
            }
        }

        struct TaskQueue {
        public:
            template<typename F>
            void Push(F&& f) {
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    queue.emplace(std::forward<F>(f));
                }
                is_ready.notify_one();
            }

            bool Pop(std::function<void()>* f) {
                std::unique_lock<std::mutex> lock(mutex);
                while (queue.empty() && !is_done) {
                    is_ready.wait(lock);
                }
                if (queue.empty()) {
                    return false;
                }
                *f = std::move(queue.front());
                queue.pop();
                return true;
            }

            template<typename F>
            bool TryPush(F&& f) {
                {
                    std::unique_lock<std::mutex> lock(mutex, std::try_to_lock);
                    if (!lock) {
                        return false;
                    }
                    queue.emplace(std::forward<F>(f));
                }
                is_ready.notify_one();
                return true;
            }

            bool TryPop(std::function<void()>* f) {
                std::unique_lock<std::mutex> lock(mutex, std::try_to_lock);
                if (!lock || queue.empty()) {
                    return false;
                }
                *f = std::move(queue.front());
                queue.pop();
                return true;
            }

            void Done() {
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    is_done = true;
                }
                is_ready.notify_all();
            }

            std::queue<std::function<void()>> queue;
            std::mutex mutex;
            std::condition_variable is_ready;
            bool is_done = false;
        };

        std::vector<std::thread> threads_;
        std::unordered_map<std::thread::id, std::size_t> thread_map_;
        std::vector<TaskQueue> queues_;
        std::atomic<std::size_t> task_id_;
    };

}  // namespace thread_pool

#endif  // THREAD_POOL_THREAD_POOL_HPP_
