#pragma once

#include <vector>
#include <memory>
#include <mutex>
#include <complex>
#include <unordered_map>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Thread-safe memory pool for grid allocations
 * 
 * This class manages a pool of pre-allocated memory buffers to avoid
 * frequent allocations/deallocations that can cause heap fragmentation
 */
class GridMemoryPool {
private:
    struct Buffer {
        std::unique_ptr<std::complex<double>[]> data;
        size_t size;
        bool in_use;
        
        Buffer(size_t sz) : size(sz), in_use(false) {
            data = std::make_unique<std::complex<double>[]>(sz);
        }
    };
    
    mutable std::mutex pool_mutex;
    std::vector<std::unique_ptr<Buffer>> buffers;
    size_t max_buffers = 10;  // Maximum number of buffers to keep
    
    // Statistics
    size_t allocations = 0;
    size_t reuses = 0;
    
public:
    /**
     * @brief Get a buffer of at least the requested size
     */
    std::complex<double>* acquire(size_t size) {
        std::lock_guard<std::mutex> lock(pool_mutex);
        
        // First, try to find an unused buffer of sufficient size
        for (auto& buf : buffers) {
            if (!buf->in_use && buf->size >= size) {
                buf->in_use = true;
                reuses++;
                // Zero out the buffer before returning
                std::fill_n(buf->data.get(), size, std::complex<double>(0.0, 0.0));
                return buf->data.get();
            }
        }
        
        // No suitable buffer found, create a new one
        if (buffers.size() < max_buffers) {
            buffers.push_back(std::make_unique<Buffer>(size));
            buffers.back()->in_use = true;
            allocations++;
            return buffers.back()->data.get();
        }
        
        // Pool is full, find the smallest buffer and resize it
        size_t smallest_idx = 0;
        size_t smallest_size = buffers[0]->size;
        for (size_t i = 1; i < buffers.size(); i++) {
            if (!buffers[i]->in_use && buffers[i]->size < smallest_size) {
                smallest_idx = i;
                smallest_size = buffers[i]->size;
            }
        }
        
        // Resize the buffer if needed
        if (buffers[smallest_idx]->size < size) {
            buffers[smallest_idx] = std::make_unique<Buffer>(size);
        }
        buffers[smallest_idx]->in_use = true;
        allocations++;
        std::fill_n(buffers[smallest_idx]->data.get(), size, std::complex<double>(0.0, 0.0));
        return buffers[smallest_idx]->data.get();
    }
    
    /**
     * @brief Release a buffer back to the pool
     */
    void release(std::complex<double>* ptr) {
        if (!ptr) return;
        
        std::lock_guard<std::mutex> lock(pool_mutex);
        
        for (auto& buf : buffers) {
            if (buf->data.get() == ptr) {
                buf->in_use = false;
                return;
            }
        }
    }
    
    /**
     * @brief Clear all buffers from the pool
     */
    void clear() {
        std::lock_guard<std::mutex> lock(pool_mutex);
        buffers.clear();
        allocations = 0;
        reuses = 0;
    }
    
    /**
     * @brief Get pool statistics
     */
    void getStats(size_t& total_allocs, size_t& total_reuses) const {
        std::lock_guard<std::mutex> lock(pool_mutex);
        total_allocs = allocations;
        total_reuses = reuses;
    }
};

// Global memory pool instance
inline GridMemoryPool& getGridMemoryPool() {
    static GridMemoryPool pool;
    return pool;
}

/**
 * @brief RAII wrapper for grid memory
 */
class GridMemoryHandle {
private:
    std::complex<double>* ptr;
    size_t size;
    
public:
    GridMemoryHandle(size_t sz) : size(sz) {
        ptr = getGridMemoryPool().acquire(sz);
    }
    
    ~GridMemoryHandle() {
        if (ptr) {
            getGridMemoryPool().release(ptr);
        }
    }
    
    // Delete copy constructor and assignment
    GridMemoryHandle(const GridMemoryHandle&) = delete;
    GridMemoryHandle& operator=(const GridMemoryHandle&) = delete;
    
    // Allow move
    GridMemoryHandle(GridMemoryHandle&& other) noexcept 
        : ptr(other.ptr), size(other.size) {
        other.ptr = nullptr;
        other.size = 0;
    }
    
    GridMemoryHandle& operator=(GridMemoryHandle&& other) noexcept {
        if (this != &other) {
            if (ptr) {
                getGridMemoryPool().release(ptr);
            }
            ptr = other.ptr;
            size = other.size;
            other.ptr = nullptr;
            other.size = 0;
        }
        return *this;
    }
    
    std::complex<double>* get() { return ptr; }
    const std::complex<double>* get() const { return ptr; }
    size_t getSize() const { return size; }
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc