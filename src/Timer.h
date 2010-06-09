#ifndef TIMER
#define TIMER
#include <ctime>
class Timer{
	clock_t startC,finishC;
	double time;
public:
	Timer(): time(0) { }
	void start(){
		startC = clock();
	}
	void stop(){
		finishC = clock();
		time = (double(finishC)-double(startC))/CLOCKS_PER_SEC;
	}
	double show() const{
		return time;
	}
	double operator()(){
		return show();
	}
};
template<class CharT, class Traits> std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const Timer& timer){

	os << timer.show();
	return os;
}
#endif 