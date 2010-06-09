#ifndef PRIORITY_QUEUE
#define PRIORITY_QUEUE

#include "heap.h"
#include <limits>

template <class T>
class PriorityQueue {
	const T MINIMO;
	Heap* heap_;
	T* cost_;
	int maxItems_;
public: 
	PriorityQueue(int maxItems, T* f) : maxItems_(maxItems), cost_(f), MINIMO(std::numeric_limits<T>::min())
	{
		heap_ = HeapInit(maxItems+1);
		HeapSetKeys(heap_, (heapvalues*)cost_);
	}
	int key_top(){
		return HeapMin(heap_);
	}
	//double value_top(){

	//}
	void pop(){
		if(!empty())
			HeapDelMin(heap_);
		else
			exit(999);
	}
	void push(int key){
		HeapInsert(heap_,key);
	}
	bool empty(){
		return HeapSize(heap_) == 0;
	}
	void decrease(int key){
		HeapDecKey(heap_, key);
	}
	void increase(int key){
		HeapIncKey(heap_, key);
	}
	void remove(int key){
		T temp = cost_[key];
		cost_[key] = MINIMO;
		decrease(key);
		
		if( key != key_top() )
			exit(100);
		pop();
		cost_[key] = temp;
	}
	void clear(){
		if(!empty())
		{
			HeapFree(heap_);
			heap_ = HeapInit(maxItems_);
			HeapSetKeys(heap_, (heapvalues*)cost_);
		}
	}
	~PriorityQueue(){
		HeapFree(heap_);
	}
};


//#include <vector>
//#include <algorithm>
//
//
//template <class T>
//class PriorityQueue {
//	const double MINIMO;
//	std::vector<int> heap_;
//	
//	int maxItems_;
//public: 
//	T* cost_;
//
//	PriorityQueue(int maxItems, T* f) : maxItems_(maxItems), cost_(f), MINIMO(std::numeric_limits<T>::min())
//	{
//
//	}
//	int key_top(){
//		T mi = 9999999;
//		int key = 0;
//		for(std::vector<int>::iterator it = heap_.begin(); it != heap_.end(); it++)
//			if(cost_[*it] < mi)
//			{
//				mi = cost_[*it];
//				key = *it;
//			}
//		return key;
//	}
//	//double value_top(){
//
//	//}
//	void pop(){
//		remove(key_top());
//	}
//	void push(int key){
//		heap_.push_back(key);
//	}
//	bool empty(){
//		return heap_.size() == 0;
//	}
//	//bool cmp(const int &a, const int &b){
//	//	//return cost_[a] >= cost_[b];
//	//	return true;
//	//}
//	void decrease(int key){
//
//	}
//	void increase(int key){
//
//	}
//	void remove(int key){
//		std::vector<int>::iterator it = std::find(heap_.begin(), heap_.end(), key);
//		if(it == heap_.end())
//			exit(key);
//		heap_.erase(it);
//	}
//	void clear(){
//		heap_.clear();
//	}
//	~PriorityQueue(){
//		heap_.clear();
//	}
//};
//
#endif 